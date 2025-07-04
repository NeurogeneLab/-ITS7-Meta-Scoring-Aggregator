# Meta Scoring Aggregator

A modular scoring engine to evaluate drug-like small molecules across multiple dimensions such as docking, ADME/PK, toxicity, druggability, drug-likeness, and off-target selectivity.

---

## 📌 Features

- 🔬 **Multi-module evaluation:** Handles six critical pharmacological modules

- ⚙️ **YAML-driven configuration:** Easily tweak weights and thresholds

- 🧪 **Comprehensive test suite:** Full coverage for edge cases, validation, config errors

- 🚨 **Override safety triggers:** Toxicity modules can auto-reject compounds

- 🔧 **Extensible scoring:** Scoring functions are modular and customizable

---

## 🏗️ Project Structure

```
meta_scorer/
├── scoring_aggregator.py       # Main scoring logic and aggregation pipeline
├── score_config.yaml           # Weights and thresholds used for scoring
├── test_scoring_aggregator.py  # Complete pytest suite
├── examples/
│   ├── input_example.json      # Input compound data
│   └── output_example.json     # Generated scoring output
```

---

## 🚀 How to Run

### 1. Install dependencies:
```bash
pip install pyyaml pytest
```

### 2. Run the scorer:
```bash
python scoring_aggregator.py
```
> Requires `examples/input_example.json` and `score_config.yaml` to exist in project.

### 3. Run unit tests:
```bash
pytest test_scoring_aggregator.py
```

---

## 📊 Example Output
```json
{
  "svs_score": 82,
  "go_decision": "go",
  "score_breakdown": {
    "docking": 20,
    "adme_pk": 20,
    "toxicity": 10,
    "druggability": 12,
    "drug_likeness": 5,
    "off_target": 15
  },
  "failure_rationale": []
}
```

---

## 🧪 Test Coverage Highlights

- ✅ High-scoring compound returns GO
- ✅ Toxicity overrides force NO-GO
- ✅ YAML config structure validation
- ✅ Missing/invalid input detection
- ✅ Scoring logic clamps outputs [0–100]

---
