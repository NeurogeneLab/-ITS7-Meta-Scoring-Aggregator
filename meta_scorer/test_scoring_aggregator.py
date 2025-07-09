import pytest
import tempfile
import os
from scoring_aggregator import aggregate_scores, MetaScoringAggregator, ScoringError


class TestMetaScoringAggregator:
    """Main test suite for MetaScoringAggregator core functionality."""

    def setup_method(self):
        """Initialize reusable valid input for tests."""
        self.valid_input = {
            "docking": {"affinity_kcal": -8.6, "rmsd": 2.1, "num_converged_poses": 4},
            "adme_pk": {
                "logP": 2.3,
                "HIA": "high",
                "half_life_hr": 3.4,
                "CYP_inhibition": False,
                "clearance": "acceptable"
            },
            "toxicity": {
                "toxicity_flag": True,
                "num_models_flagged": 1,
                "high_risk_models": []
            },
            "druggability": {
                "pocket_score": 0.78,
                "volume": 420,
                "hydrophobicity": 0.6,
                "polarity": 0.3
            },
            "drug_likeness": {"num_violations": 1, "passes": True},
            "off_target": {
                "num_off_targets": 4,
                "avg_affinity_offtarget": -6.1,
                "selectivity_index": 6.5
            }
        }

    def test_high_scoring_go(self):
        """Should return 'go' for near-perfect compound with max scores."""
        input_data = {
            "docking": {"affinity_kcal": -9.5, "rmsd": 1.5, "num_converged_poses": 5},
            "adme_pk": {
                "logP": 3.0,
                "HIA": "high",
                "half_life_hr": 5,
                "CYP_inhibition": False,
                "clearance": "acceptable"
            },
            "toxicity": {
                "toxicity_flag": False,
                "num_models_flagged": 0,
                "high_risk_models": []
            },
            "druggability": {
                "pocket_score": 0.8,
                "volume": 450,
                "hydrophobicity": 0.6,
                "polarity": 0.2
            },
            "drug_likeness": {"num_violations": 0, "passes": True},
            "off_target": {
                "num_off_targets": 3,
                "avg_affinity_offtarget": -5.0,
                "selectivity_index": 12
            }
        }

        result = aggregate_scores(input_data)

        # All scores max out, total 100 (clamped)
        assert result["go_decision"] == "go"
        assert result["svs_score"] == 100
        assert result["failure_rationale"] == []

        # Confirm each module's score
        expected_breakdown = {
            "docking": 25,
            "adme_pk": 20,
            "toxicity": 20,
            "druggability": 12,
            "drug_likeness": 10,
            "off_target": 18
        }
        assert result["score_breakdown"] == expected_breakdown

    def test_borderline_still_go(self):
        """Test borderline compound just clearing 70+ mark."""
        input_data = self.valid_input.copy()
        input_data["drug_likeness"] = {"num_violations": 1, "passes": True}
        result = aggregate_scores(input_data)
        assert result["go_decision"] == "go"
        assert result["svs_score"] >= 70

    def test_toxic_no_go_override(self):
        """Should return 'no-go' if toxicity override is triggered."""
        input_data = self.valid_input.copy()
        input_data["toxicity"] = {
            "toxicity_flag": True,
            "num_models_flagged": 4,  # triggers override
            "high_risk_models": ["DILIrank", "ProTox-II"]
        }
        result = aggregate_scores(input_data)
        assert result["go_decision"] == "no-go"
        assert any("toxicity" in reason.lower() for reason in result["failure_rationale"])

    def test_low_score_no_go(self):
        """Should return 'no-go' if total score is below threshold."""
        input_data = {
            "docking": {"affinity_kcal": -5.0, "rmsd": 4.0, "num_converged_poses": 1},
            "adme_pk": {
                "logP": 8.0,
                "HIA": "low",
                "half_life_hr": 1.0,
                "CYP_inhibition": True,
                "clearance": "poor"
            },
            "toxicity": {
                "toxicity_flag": True,
                "num_models_flagged": 2,
                "high_risk_models": ["ProTox-II"]
            },
            "druggability": {
                "pocket_score": 0.3,
                "volume": 200,
                "hydrophobicity": 0.2,
                "polarity": 0.8
            },
            "drug_likeness": {"num_violations": 3, "passes": False},
            "off_target": {
                "num_off_targets": 20,
                "avg_affinity_offtarget": -7.0,
                "selectivity_index": 1.5
            }
        }
        result = aggregate_scores(input_data)
        assert result["go_decision"] == "no-go"
        assert result["svs_score"] < 70

    def test_high_risk_toxicity_models(self):
        """Triggers penalty from known high-risk models even if model count is low."""
        input_data = self.valid_input.copy()
        input_data["toxicity"] = {
            "toxicity_flag": True,
            "num_models_flagged": 1,
            "high_risk_models": ["DILIrank", "ProTox-II"]
        }
        result = aggregate_scores(input_data)
        assert result["score_breakdown"]["toxicity"] <= 0
        assert result["go_decision"] == "no-go"
        
    def test_edge_case_boundary_values(self):
        """Test boundary values for scoring thresholds (inclusive logic)."""
        input_data = {
              "docking": {"affinity_kcal": -9.0, "rmsd": 3.0, "num_converged_poses": 3},
              "adme_pk": {
                  "logP": 0,
                  "HIA": "HIGH",
                  "half_life_hr": 2.0,
                  "CYP_inhibition": False,
                  "clearance": "ACCEPTABLE"
              },
              "toxicity": {
                  "toxicity_flag": True,
                  "num_models_flagged": 1,  # ✅ Avoids override
                  "high_risk_models": []
              },
              "druggability": {
                  "pocket_score": 0.7,  # ✅ Right at high threshold
                  "volume": 400,
                  "hydrophobicity": 0.5,
                  "polarity": 0.3
              },
              "drug_likeness": {"num_violations": 2, "passes": False},
              "off_target": {
                  "num_off_targets": 5,
                  "avg_affinity_offtarget": -6.0,
                  "selectivity_index": 10
              }
          }

        result = aggregate_scores(input_data)
        assert result["go_decision"] == "go"
        assert result["svs_score"] == 70


    def test_missing_optional_fields(self):
        """Ensure defaults handle missing optional fields gracefully."""
        input_data = {
            "docking": {"affinity_kcal": -8.0},
            "adme_pk": {"logP": 2.0},
            "toxicity": {},
            "druggability": {"pocket_score": 0.8},
            "drug_likeness": {},
            "off_target": {"selectivity_index": 5.0}
        }
        result = aggregate_scores(input_data)
        assert "svs_score" in result
        assert result["go_decision"] in ["go", "no-go"]


class TestInputValidation:
    """Test robustness of input validation."""

    def test_invalid_input_type(self):
        with pytest.raises(ScoringError, match="Input must be a dictionary"):
            aggregate_scores("invalid")

    def test_empty_input(self):
        with pytest.raises(ScoringError, match="Input dictionary cannot be empty"):
            aggregate_scores({})

    def test_missing_required_modules(self):
        incomplete_input = {
            "docking": {"affinity_kcal": -8.0},
            "adme_pk": {"logP": 2.0}
        }
        with pytest.raises(ScoringError, match="Missing required modules"):
            aggregate_scores(incomplete_input)


class TestConfigHandling:
    """Ensure YAML config loading and structure validation works."""

    def test_missing_config_file(self):
        with pytest.raises(ScoringError, match="Configuration file not found"):
            MetaScoringAggregator("fake_config.yaml")

    def test_invalid_yaml_config(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write("invalid: yaml: [")
            bad_path = f.name
        try:
            with pytest.raises(ScoringError, match="Invalid YAML format"):
                MetaScoringAggregator(bad_path)
        finally:
            os.unlink(bad_path)

    def test_empty_config_file(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write("")
            bad_path = f.name
        try:
            with pytest.raises(ScoringError, match="Configuration file is empty"):
                MetaScoringAggregator(bad_path)
        finally:
            os.unlink(bad_path)

    def test_incomplete_config_structure(self):
        bad_config = """
        weights:
          docking: 25
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            f.write(bad_config)
            path = f.name
        try:
            with pytest.raises(ScoringError, match="Missing required config section"):
                MetaScoringAggregator(path)
        finally:
            os.unlink(path)


class TestScoreCalculationAccuracy:
    """Confirm precise breakdown logic matches expectations."""
        
    def test_example_calculation_accuracy(self):
        input_data = {
            "docking": {"affinity_kcal": -8.6, "rmsd": 2.1, "num_converged_poses": 4},
            "adme_pk": {
                "logP": 2.3,
                "HIA": "high",
                "half_life_hr": 3.4,
                "CYP_inhibition": False,
                "clearance": "acceptable"
            },
            "toxicity": {
                "toxicity_flag": True,
                "num_models_flagged": 1,
                "high_risk_models": []
            },
            "druggability": {
                "pocket_score": 0.78,
                "volume": 420,
                "hydrophobicity": 0.6,
                "polarity": 0.3
            },
            "drug_likeness": {"num_violations": 1, "passes": True},
            "off_target": {
                "num_off_targets": 4,
                "avg_affinity_offtarget": -6.1,
                "selectivity_index": 6.5
            }
        }

        result = aggregate_scores(input_data)

        expected_scores = {
            "docking": 20,
            "adme_pk": 20,
            "toxicity": 10,
            "druggability": 12,
            "drug_likeness": 5,
            "off_target": 13
        }

        assert result["score_breakdown"] == expected_scores
        assert result["svs_score"] == 80
        assert result["go_decision"] == "go"


    def test_score_clamping(self):
        """Scores must be clamped between 0–100."""
        bad_input = {
            "docking": {"affinity_kcal": 0, "rmsd": 10, "num_converged_poses": 0},
            "adme_pk": {
                "logP": 10,
                "HIA": "low",
                "half_life_hr": 0,
                "CYP_inhibition": True,
                "clearance": "poor"
            },
            "toxicity": {
                "num_models_flagged": 2,
                "high_risk_models": ["DILIrank", "ProTox-II"]
            },
            "druggability": {"pocket_score": 0.1, "volume": 100, "hydrophobicity": 0.1},
            "drug_likeness": {"num_violations": 5, "passes": False},
            "off_target": {"num_off_targets": 50, "selectivity_index": 0.5}
        }
        result = aggregate_scores(bad_input)
        assert 0 <= result["svs_score"] <= 100
        
    def test_score_just_below_go(self):
        """Test a compound with total score 69 and no override – should be no-go."""
        input_data = {
            "docking": {"affinity_kcal": -7.6, "rmsd": 2.9, "num_converged_poses": 3},
            "adme_pk": {
                "logP": 5.0,
                "HIA": "high",
                "half_life_hr": 2.0,
                "CYP_inhibition": False,
                "clearance": "acceptable"
            },
            "toxicity": {
                "toxicity_flag": True,
                "num_models_flagged": 1,
                "high_risk_models": []
            },
            "druggability": {
                "pocket_score": 0.6,
                "volume": 390,
                "hydrophobicity": 0.4,
                "polarity": 0.3
            },
            "drug_likeness": {"num_violations": 1, "passes": True},
            "off_target": {
                "num_off_targets": 5,
                "avg_affinity_offtarget": -6.0,
                "selectivity_index": 4.9
            }
        }
        result = aggregate_scores(input_data)
        assert result["go_decision"] == "no-go"
        assert result["svs_score"] == 69

    def test_exact_thresholds_inclusive(self):
        """Test that exact threshold values receive proper scores with inclusive logic."""
        input_data = {
            "docking": {"affinity_kcal": -9.0, "rmsd": 3.0, "num_converged_poses": 3},
            "adme_pk": {
                "logP": 0,
                "HIA": "HIGH",
                "half_life_hr": 2.0,
                "CYP_inhibition": False,
                "clearance": "ACCEPTABLE"
            },
            "toxicity": {
                "toxicity_flag": True,
                "num_models_flagged": 1,
                "high_risk_models": []
            },
            "druggability": {
                "pocket_score": 0.7,
                "volume": 400,
                "hydrophobicity": 0.5,
                "polarity": 0.3
            },
            "drug_likeness": {"num_violations": 2, "passes": False},
            "off_target": {
                "num_off_targets": 5,
                "avg_affinity_offtarget": -6.0,
                "selectivity_index": 10
            }
        }
        result = aggregate_scores(input_data)
        assert result["score_breakdown"]["docking"] == 25
        assert result["score_breakdown"]["druggability"] <= 10
        assert result["go_decision"] == "go"

    def test_module_score_clamping(self):
        """Ensure that individual modules don't exceed their max weights."""
        input_data = {
            "docking": {"affinity_kcal": -10.0, "rmsd": 1.0, "num_converged_poses": 5},
            "adme_pk": {
                "logP": 2.5,
                "HIA": "high",
                "half_life_hr": 6.0,
                "CYP_inhibition": False,
                "clearance": "acceptable"
            },
            "toxicity": {
                "toxicity_flag": False,
                "num_models_flagged": 0,
                "high_risk_models": []
            },
            "druggability": {
                "pocket_score": 0.9,
                "volume": 500,
                "hydrophobicity": 0.6,
                "polarity": 0.3
            },
            "drug_likeness": {"num_violations": 0, "passes": True},
            "off_target": {
                "num_off_targets": 2,
                "avg_affinity_offtarget": -5.5,
                "selectivity_index": 20
            }
        }
        result = aggregate_scores(input_data)
        assert result["score_breakdown"]["off_target"] <= 15  # weight limit
        assert result["score_breakdown"]["druggability"] <= 10  # weight limit

