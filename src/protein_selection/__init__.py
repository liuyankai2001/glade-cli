"""Local protein sequence acquisition, independent of enzyme selection."""

from .download import download_protein_sequence
from .add import add_protein_to_reaction

__all__ = ["download_protein_sequence", "add_protein_to_reaction"]
