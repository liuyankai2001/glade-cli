"""Local protein sequence acquisition, independent of enzyme selection."""

from .download import download_protein_sequence
from .add import add_protein_to_reaction
from .remove import remove_protein_from_reaction
from .replace import replace_protein_in_reaction

__all__ = ["download_protein_sequence", "add_protein_to_reaction", "remove_protein_from_reaction", "replace_protein_in_reaction"]
