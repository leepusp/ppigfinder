"""
AlphaFold 3 JSON Builder Module.
Generates compliant JSON input files for AlphaFold 3 structural predictions.
"""
import json
import logging
from typing import List, Dict, Any

logger = logging.getLogger(__name__)

class AF3JobBuilder:
    def __init__(self, job_name: str, seed: int = 1):
        """Initialises the AlphaFold 3 job builder."""
        self.job_name = job_name
        self.seed = seed
        self.sequences = []

    def add_protein_sequence(self, sequence: str, copies: int = 1) -> None:
        """Adds a protein sequence to the complex."""
        self.sequences.append({
            "protein": {
                "id": [str(i) for i in range(len(self.sequences) + 1, len(self.sequences) + copies + 1)],
                "sequence": sequence.upper()
            }
        })

    def build_json_dict(self) -> List[Dict[str, Any]]:
        """Constructs the final dictionary structure for AF3."""
        if not self.sequences:
            logger.warning("No sequences added to the AF3 job.")
            
        return [{
            "name": self.job_name,
            "modelSeeds": [self.seed],
            "sequences": self.sequences,
            "dialect": "alphafold3",
            "version": 1
        }]

    def save_json(self, filepath: str) -> None:
        """Exports the job configuration to a JSON file."""
        job_data = self.build_json_dict()
        try:
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(job_data, f, indent=4)
            logger.info(f"AF3 JSON saved to {filepath}")
        except IOError as e:
            logger.error(f"Failed to save AF3 JSON to {filepath}: {e}")
            raise
