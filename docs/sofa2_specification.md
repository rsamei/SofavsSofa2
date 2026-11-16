# SOFA-2 implementation in MIMIC-IV v3.1

Each component against the published definition (Ranzani et al., JAMA 2025, Table 2 and footnotes). Scores are computed
for every ICU hour (`sql/scores/01_sofa2_withintime.sql`, `02_cv_hourly.sql`, `03_sofa2_hourly.sql`). The score at an
hour is the worst value of each component over the preceding 24 hours. A component with no measurement scores 0.
Concept tables are mimic-code at commit 4582f39 (`sql/mimic_code/`).

## Respiratory

| Points | SOFA-2 | Implementation |
|---|---|---|
| 0 | PaO2/FiO2 > 300 | as published |
| 1 | ≤ 300 | as published |
| 2 | ≤ 225 | as published |
| 3 | ≤ 150 with advanced respiratory support | as published |
| 4 | ≤ 75 with advanced respiratory support, or ECMO | as published |
| SpO2/FiO2 (footnote) | used when PaO2/FiO2 is unavailable and SpO2 < 98%: > 300, ≤ 300, ≤ 250, ≤ 200 and ≤ 120 with support | as published |

- **PaO2/FiO2:** the lowest value in the hour from arterial samples in the `bg` concept (lab items 50821 PO2 and 50816 FiO2). When the blood gas has no FiO2, `bg` takes FiO2 charted up to 4 hours before (chart item 223835).
- **SpO2/FiO2:**
  - SpO2 from `vitalsign` (item 220277), only values below 98%.
  - FiO2 from `ventilator_setting` (item 223835), accepted only between 0.21 and 1.0, and matched only when charted at or up to 60 minutes before the SpO2.
  - FiO2 is charted during ventilation and high-flow oxygen, and also with oxygen masks and other devices for which a set or estimated fraction is recorded (for example aerosol mask, face tent, Venturi mask, non-rebreather and, less often, nasal cannula); hours without a charted FiO2, including most room-air hours, have no ratio.
- **Advanced respiratory support:** `ventilation` status InvasiveVent, NonInvasiveVent, HFNC or Tracheostomy.
  - For PaO2/FiO2 it is taken at the blood-gas time.
  - For SpO2/FiO2 it is any support episode overlapping the hour.
- **ECMO:** chart item 229268 (circuit configuration VV, VA or VAV) in the hour.

## Brain

| Points | SOFA-2 | Implementation |
|---|---|---|
| 0–4 | GCS 15 / 13–14 / 9–12 / 6–8 / 3–5 | as published, on the lowest total GCS of the hour |
| Footnote d | if not all three domains can be assessed, use the best motor response: 6 / 5 / 4 / 3 / ≤ 2 | verbal "No Response-ETT" marks the total as not assessable; the best motor score of the hour is used |
| Footnote c | sedated patients: use the last GCS before sedation | sedation = propofol, midazolam, lorazepam or dexmedetomidine in `icu.inputevents`. Hours within 6 hours on both sides of a sedated hour are bridged into one sedation block. The last assessable pre-sedation score is carried forward for up to 72 hours; with none, the score is 0 |
| Footnote e | at least 1 point when receiving drug therapy for delirium | haloperidol, quetiapine, olanzapine, risperidone, ziprasidone or chlorpromazine administered in `hosp.emar` during the ICU stay, in the hour |

GCS comes from chart items 220739 (eyes), 223900 (verbal) and 223901 (motor).

## Cardiovascular

| Points | SOFA-2 | Implementation |
|---|---|---|
| 0 | MAP ≥ 70 mmHg, no vasopressor or inotrope | as published |
| 1 | MAP < 70 mmHg, no vasopressor or inotrope | as published (lowest MAP in the hour) |
| 2 | norepinephrine + epinephrine ≤ 0.2 µg/kg/min, or any dose of another vasopressor or inotrope | as published |
| 3 | norepinephrine + epinephrine > 0.2–0.4, or ≤ 0.2 with another vasopressor or inotrope | as published |
| 4 | norepinephrine + epinephrine > 0.4, or > 0.2–0.4 with another agent, or mechanical circulatory support | as published |
| Dopamine alone | ≤ 20 / > 20–40 / > 40 µg/kg/min: 2 / 3 / 4 | as published |

- **Doses:** infusion rates from the mimic-code medication concepts (items 221906 norepinephrine, 221289 epinephrine, 221662 dopamine, 221653 dobutamine, 221749 phenylephrine, 222315 vasopressin, 221986 milrinone), in µg/kg/min as charted.
  - MIMIC-IV records norepinephrine as the base dose used in US practice, and no salt conversion is applied.
  - The rate in an hour is the highest rate of any infusion row overlapping the hour.
- **Other vasopressor or inotrope:** vasopressin, phenylephrine, dobutamine, milrinone, or dopamine given together with norepinephrine or epinephrine.
- **Minimum duration:** an infusion counts only if its continuous episode lasts at least 60 minutes, rows being merged when the next starts within 1 minute. At the prediction horizon it counts only once 60 minutes have elapsed.
- **Mechanical circulatory support:**
  - VA or VAV ECMO (item 229268);
  - any chart item whose label contains IABP, Impella, LVAD or " VAD".

MAP is from `vitalsign` (items 220052, 220181, 225312).

## Liver

| Points | SOFA-2 | Implementation |
|---|---|---|
| 0–4 | bilirubin ≤ 1.2 / ≤ 3.0 / ≤ 6.0 / ≤ 12.0 / > 12.0 mg/dL | as published, highest total bilirubin of the hour (`enzyme`, lab item 50885) |

## Hemostasis

| Points | SOFA-2 | Implementation |
|---|---|---|
| 0–4 | platelets > 150 / ≤ 150 / ≤ 100 / ≤ 80 / ≤ 50 × 10³/µL | as published, lowest platelet count of the hour (`complete_blood_count`, lab item 51265) |

## Kidney

| Points | SOFA-2 | Implementation |
|---|---|---|
| 0 | creatinine ≤ 1.20 mg/dL | as published (`chemistry`, lab item 50912) |
| 1 | creatinine ≤ 2.0, or urine output < 0.5 mL/kg/h for 6–12 h | see urine output below |
| 2 | creatinine ≤ 3.5, or urine output < 0.5 mL/kg/h for ≥ 12 h | as above |
| 3 | creatinine > 3.5, or urine output < 0.3 mL/kg/h for ≥ 24 h, or anuria for ≥ 12 h | as above |
| 4 | receiving or fulfilling criteria for RRT | RRT from the `rrt` concept. Criteria: creatinine > 1.2 mg/dL or urine output < 0.3 mL/kg/h for 6 h, together with potassium ≥ 6.0 mmol/L or metabolic acidosis (pH ≤ 7.20 and bicarbonate ≤ 12 mmol/L) |

- **Urine output, primary rule:**
  - Source: the `urine_output` concept (items 226557–226567, 226584, 227488, 227489).
  - Each chart's collection interval is the time since the previous urine chart; for the first chart it is the time since ICU admission.
  - The rate over a 6-, 12- or 24-hour window is the volume charted in the window divided by the summed collection intervals of those charts and by weight, as in mimic-code `urine_output_rate`.
  - A window is scored only when those intervals cover the whole window.
  - Anuria needs zero volume over a covered 12-hour window with at least two charts.
- **Weight:**
  - Source: chart items 226512 (admission weight) and 224639 (daily weight).
  - The weight at an hour is the last weight charted by then.
  - Before the first charted weight, the first admission weight is back-applied only if it was charted by the end of the prediction horizon.
- **Variant s (single-chart rule):** the volume charted in the trailing window divided by the window length and weight, so one chart can score a window.
- **Variant g (coverage rule):** the same as variant s, but with at least two charts and charting intervals covering 5 of 6, 10 of 12 or 22 of 24 hours.
- **Variant w (bounded weight):** the primary rule, but hours whose weight is below 30 kg or above 300 kg lose the urine criteria (creatinine and RRT only).

## Time windows

- The day-1 score is the 24-hour component-worst total at the end of ICU hour 23, or at the last scored hour if the stay ends earlier.
- Laboratory values are matched to hours by hospital admission, so for stays shorter than 24 hours the day-1 window can include pre-ICU hours for four components.
- Look-ahead is truncated at the horizon for sedation bridging, back-applied weight, the SpO2/FiO2 FiO2 match and the 60-minute infusion rule.
