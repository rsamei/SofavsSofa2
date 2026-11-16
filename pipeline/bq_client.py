"""Minimal BigQuery REST client and configuration for the pipeline (Python standard library only).

Configuration (environment variables):
  BQ_PROJECT        Google Cloud project that runs the queries and holds the datasets (required).
  UPSTREAM_DATASET  dataset for the mimic-code concepts, cohort_base and sofa_withintime (default sofa2_upstream).
  BUILD_DATASET     dataset for the hourly scores and the cohort (default sofa2_build).
  OUTPUT_DIR        folder for the JSON results (default output).

Credentials, in this order: GCP_ACCESS_TOKEN or CLOUDSDK_AUTH_ACCESS_TOKEN (e.g. `gcloud auth print-access-token`),
or an authorized_user application-default credential (`gcloud auth application-default login`). The account needs
PhysioNet access to physionet-data (MIMIC-IV v3.1).

SQL may use the placeholders @UPSTREAM@ and @BUILD@, replaced by project.dataset before each query.
"""
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
from pathlib import Path

API = "https://bigquery.googleapis.com/bigquery/v2"
TOKEN_URL = "https://oauth2.googleapis.com/token"

PROJECT = os.environ.get("BQ_PROJECT", "")
UPSTREAM_DATASET = os.environ.get("UPSTREAM_DATASET", "sofa2_upstream")
BUILD_DATASET = os.environ.get("BUILD_DATASET", "sofa2_build")
UPSTREAM, BUILD = "@UPSTREAM@", "@BUILD@"
OUT = Path(os.environ.get("OUTPUT_DIR", "output"))
OUT.mkdir(parents=True, exist_ok=True)
ROOT = Path(__file__).resolve().parent.parent
SQL_DIR = ROOT / "sql"


class BQError(RuntimeError):
    pass


def resolve(sql):
    if not PROJECT:
        raise BQError("set BQ_PROJECT")
    return sql.replace(UPSTREAM, f"{PROJECT}.{UPSTREAM_DATASET}").replace(BUILD, f"{PROJECT}.{BUILD_DATASET}")


def _adc_info():
    path = os.environ.get("GOOGLE_APPLICATION_CREDENTIALS") or str(
        Path.home() / ".config/gcloud/application_default_credentials.json")
    if Path(path).is_file():
        return json.loads(Path(path).read_text(encoding="utf-8"))
    return None


_cached = {"token": None, "at": 0.0}


def access_token():
    for var in ("GCP_ACCESS_TOKEN", "CLOUDSDK_AUTH_ACCESS_TOKEN"):
        if os.environ.get(var):
            return os.environ[var].strip()
    info = _adc_info()
    if not info:
        raise BQError("no BigQuery credential: set GCP_ACCESS_TOKEN or run `gcloud auth application-default login`")
    if info.get("type") != "authorized_user":
        raise BQError(f"unsupported credential type {info.get('type')!r}; use an authorized_user credential")
    if _cached["token"] and time.time() - _cached["at"] < 2700:
        return _cached["token"]
    body = urllib.parse.urlencode({
        "client_id": info["client_id"], "client_secret": info["client_secret"],
        "refresh_token": info["refresh_token"], "grant_type": "refresh_token"}).encode()
    with urllib.request.urlopen(urllib.request.Request(TOKEN_URL, data=body), timeout=60) as r:
        _cached.update(token=json.load(r)["access_token"], at=time.time())
    return _cached["token"]


def _call(method, url, token, payload=None, retries=4):
    data = json.dumps(payload).encode() if payload is not None else None
    for attempt in range(retries + 1):
        req = urllib.request.Request(url, data=data, method=method, headers={
            "Authorization": "Bearer " + token, "Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            detail = e.read().decode("utf-8", "replace")[:500]
            # 502/503 are transient; the caller-supplied jobId makes a retried submission idempotent.
            if e.code in (502, 503) and attempt < retries:
                time.sleep(5 * (attempt + 1))
                continue
            raise BQError(f"BigQuery HTTP {e.code}: {detail}") from None


def _run_job(project, sql, dry=False, max_wait_s=3600):
    """Submit a query job (jobs.insert) and poll until it finishes."""
    token = access_token()
    job_id = "sofa2_" + uuid.uuid4().hex
    body = {"jobReference": {"projectId": project, "jobId": job_id, "location": "US"},
            "configuration": {"dryRun": dry, "query": {"query": resolve(sql), "useLegacySql": False}}}
    try:
        job = _call("POST", f"{API}/projects/{project}/jobs", token, body)
    except BQError as e:
        if "HTTP 409" not in str(e):  # 409: the job was already created by an earlier attempt
            raise
        job = _call("GET", f"{API}/projects/{project}/jobs/{job_id}?location=US", token)
    if dry:
        return job
    waited = 0
    while job.get("status", {}).get("state") != "DONE":
        if waited > max_wait_s:
            raise BQError(f"job {job_id} did not finish in {max_wait_s} s")
        time.sleep(10)
        waited += 10
        job = _call("GET", f"{API}/projects/{project}/jobs/{job_id}?location=US", access_token())
    err = job.get("status", {}).get("errorResult")
    if err:
        raise BQError(f"job {job_id} failed: {err}")
    return job


def ensure_dataset(project, dataset, location="US"):
    """Create project.dataset if it does not exist."""
    token = access_token()
    try:
        _call("GET", f"{API}/projects/{project}/datasets/{dataset}", token)
        return False
    except BQError as e:
        if "HTTP 404" not in str(e):
            raise
    _call("POST", f"{API}/projects/{project}/datasets", token,
          {"datasetReference": {"projectId": project, "datasetId": dataset}, "location": location})
    return True


def dry_run(project, sql):
    """Validate a statement without running it; returns the bytes it would process."""
    job = _run_job(project, sql, dry=True)
    return int(job.get("statistics", {}).get("totalBytesProcessed", 0))


def execute(project, sql, max_wait_s=3600):
    """Run a DDL/DML statement (e.g. CREATE TABLE AS) and return (job_id, affected_rows)."""
    job = _run_job(project, sql, max_wait_s=max_wait_s)
    q = job.get("statistics", {}).get("query", {})
    return job["jobReference"]["jobId"], q.get("numDmlAffectedRows")


def query(project, sql, max_wait_s=1800):
    """Run a GoogleSQL query and return (field_names, rows, job_id). Aggregates only."""
    job = _run_job(project, sql, max_wait_s=max_wait_s)
    job_id = job["jobReference"]["jobId"]
    res = _call("GET", f"{API}/projects/{project}/queries/{job_id}?location=US&maxResults=1000", access_token())
    fields = [f["name"] for f in res["schema"]["fields"]]
    rows = [[c.get("v") for c in r["f"]] for r in res.get("rows", [])]
    if res.get("pageToken"):
        raise BQError("result paginated; aggregate queries must return one page")
    return fields, rows, job_id
