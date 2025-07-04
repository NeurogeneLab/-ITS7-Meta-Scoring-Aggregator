import yaml
import os
import logging
from typing import Dict, List, Union, Any

# --- Logging Setup ---
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())  # Prevent logging unless explicitly configured


class ScoringError(Exception):
    """Custom exception for scoring-related errors."""
    pass


class MetaScoringAggregator:
    """Primary class responsible for aggregating scores across evaluation modules."""

    def __init__(self, config_path: str = "score_config.yaml"):
        """Initialize and load configuration file."""
        logger.info(f"Initializing MetaScoringAggregator with config: {config_path}")
        self.config = self._load_config(config_path)
        self._validate_config()

    def _load_config(self, config_path: str) -> Dict:
        """Load YAML configuration and return as a dictionary."""
        try:
            if not os.path.exists(config_path):
                logger.error(f"Configuration file not found: {config_path}")
                raise FileNotFoundError(f"Configuration file not found: {config_path}")

            with open(config_path, "r") as file:
                config = yaml.safe_load(file)

            if not config:
                raise ScoringError("Configuration file is empty or invalid")

            logger.info("Configuration loaded successfully.")
            return config

        except yaml.YAMLError as e:
            raise ScoringError(f"Invalid YAML format: {e}")
        except Exception as e:
            raise ScoringError(f"Error loading configuration: {e}")

    def _validate_config(self):
        """Ensure required sections and modules exist in the configuration."""
        required_sections = ["weights", "thresholds"]
        required_modules = ["docking", "adme_pk", "toxicity", "druggability", "drug_likeness", "off_target"]

        for section in required_sections:
            if section not in self.config:
                raise ScoringError(f"Missing required config section: {section}")

        for module in required_modules:
            if module not in self.config["weights"]:
                raise ScoringError(f"Missing weight for module: {module}")
            if module not in self.config["thresholds"]:
                raise ScoringError(f"Missing thresholds for module: {module}")

    def aggregate_scores(self, module_outputs: Dict[str, Any]) -> Dict[str, Union[int, str, Dict, List]]:
        """Main scoring function: aggregates scores from all modules."""
        logger.info("Aggregating scores for compound.")
        self._validate_input(module_outputs)

        thresholds = self.config["thresholds"]
        score_breakdown = {}
        failure_rationale = []
        total_score = 0

        # Scoring function dispatcher
        scoring_functions = {
            "docking": self._score_docking,
            "adme_pk": self._score_adme_pk,
            "toxicity": self._score_toxicity,
            "druggability": self._score_druggability,
            "drug_likeness": self._score_drug_likeness,
            "off_target": self._score_off_target,
        }

        for module_name, scoring_func in scoring_functions.items():
            module_data = module_outputs.get(module_name, {})
            module_config = thresholds[module_name]

            if module_name == "toxicity":
                score, failure = scoring_func(module_data, module_config)
                if failure:
                    logger.warning(f"Toxicity override triggered: {failure}")
                    failure_rationale.append(f"Toxicity override triggered: {failure}")
            else:
                score = scoring_func(module_data, module_config)

            score_breakdown[module_name] = score
            total_score += score

        svs_score = min(max(int(total_score), 0), 100)  # Clamp between 0–100
        go_decision = "go" if svs_score >= 70 and not failure_rationale else "no-go"

        return {
            "svs_score": svs_score,
            "go_decision": go_decision,
            "score_breakdown": score_breakdown,
            "failure_rationale": failure_rationale,
        }

    def _validate_input(self, module_outputs: Dict[str, Any]):
        """Validate required modules are present in input."""
        if not isinstance(module_outputs, dict):
            raise ScoringError("Input must be a dictionary")
        if not module_outputs:
            raise ScoringError("Input dictionary cannot be empty")

        required_modules = ["docking", "adme_pk", "toxicity", "druggability", "drug_likeness", "off_target"]
        missing = [m for m in required_modules if m not in module_outputs]
        if missing:
            raise ScoringError(f"Missing required modules: {missing}")

    # --- Module Scoring Functions ---

    def _score_docking(self, docking_data: Dict, config: Dict) -> int:
        affinity = docking_data.get("affinity_kcal", 0)
        rmsd = docking_data.get("rmsd", float('inf'))
        poses = docking_data.get("num_converged_poses", 0)
        score = 0

        if affinity < config["affinity"]["high"]:
            score = 25
        elif affinity < config["affinity"]["medium"]:
            score = 20
        elif affinity < config["affinity"]["low"]:
            score = 10

        if rmsd > config["rmsd_max"] or poses < config["min_poses"]:
            score -= 5

        return max(score, 0)

    def _score_adme_pk(self, adme_data: Dict, config: Dict) -> int:
        score = 0
        if config["logP_range"][0] <= adme_data.get("logP", -1) <= config["logP_range"][1]:
            score += 5
        if adme_data.get("HIA", "").lower() == "high":
            score += 5
        if adme_data.get("half_life_hr", 0) > config["half_life_min"]:
            score += 5
        if adme_data.get("clearance", "").lower() in [c.lower() for c in config["clearance_acceptable"]]:
            score += 5
        if adme_data.get("CYP_inhibition", False):
            score += config["cyp_inhibition_penalty"]
        return max(score, 0)

    def _score_toxicity(self, tox_data: Dict, config: Dict) -> tuple:
        score, failure_reason = 0, None
        num_flagged = tox_data.get("num_models_flagged", 0)

        if num_flagged == 0:
            score = 20
        elif num_flagged == 1:
            score = 10
        elif num_flagged > config["max_flags_allowed"]:
            failure_reason = f">{config['max_flags_allowed']} toxicity models flagged"

        high_risk_penalties = [m for m in tox_data.get("high_risk_models", []) if m in config["override_models"]]
        if high_risk_penalties:
            score -= 10 * len(high_risk_penalties)
            penalty_reason = f"High-risk models flagged: {', '.join(high_risk_penalties)}"
            failure_reason = f"{failure_reason}; {penalty_reason}" if failure_reason else penalty_reason

        return max(score, 0), failure_reason

    def _score_druggability(self, drugg_data: Dict, config: Dict) -> int:
        score = 0
        if drugg_data.get("pocket_score", 0) > config["high_score"]:
            score += 10
        elif drugg_data.get("pocket_score", 0) >= config["mid_score"]:
            score += 5

        if drugg_data.get("volume", 0) >= config["bonus_volume_min"] and \
           drugg_data.get("hydrophobicity", 0) >= config["bonus_hydrophobicity_min"]:
            score += config["volume_hydrophobicity_bonus"]

        return score

    def _score_drug_likeness(self, dl_data: Dict, config: Dict) -> int:
        violations = dl_data.get("num_violations", float('inf'))
        score_table = config["score_by_violations"]
        violations = min(violations, max(score_table.keys()))
        return score_table.get(violations, 0)

    def _score_off_target(self, ot_data: Dict, config: Dict) -> int:
        si = ot_data.get("selectivity_index", 0)
        thresholds = config["si_thresholds"]
        score = 0

        if si > thresholds[2]:
            score = 15
        elif si > thresholds[1]:
            score = 10
        elif si > thresholds[0]:
            score = 5

        if ot_data.get("num_off_targets", float('inf')) < config["low_target_bonus_threshold"]:
            score += config["low_target_bonus"]

        return score


# Public interface function with logging setup

def aggregate_scores(module_outputs: Dict[str, Any]) -> Dict[str, Union[int, str, Dict, List]]:
    """Run scoring aggregation with default config and logging."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )
    aggregator = MetaScoringAggregator()
    return aggregator.aggregate_scores(module_outputs)


# --- Example CLI usage ---

if __name__ == "__main__":
    import json

    input_path = "examples/input_example.json"
    output_path = "examples/output_example.json"

    if not os.path.exists(input_path):
        print(f"❌ Cannot find input file: {input_path}")
    else:
        with open(input_path) as f:
            input_data = json.load(f)

        result = aggregate_scores(input_data)

        print("\n✅ Scoring result:")
        print(json.dumps(result, indent=2))

        os.makedirs("examples", exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(result, f, indent=2)

        print(f"\n📁 Output saved to: {output_path}")
