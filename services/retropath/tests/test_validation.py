import pytest

from app.models import JobParameters
from app.validation import InputValidationError, validate_compound_csv


SOURCE = b'Name,InChI\n"target","InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)"\n'
SINK = b'Name,InChI\n"C00022","InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)"\n'


def test_valid_source_and_sink():
    assert validate_compound_csv(SOURCE, kind="source").row_count == 1
    assert validate_compound_csv(SINK, kind="sink").row_count == 1


@pytest.mark.parametrize(
    "payload,kind",
    [
        (b"wrong,header\na,b\n", "source"),
        (b"Name,InChI\n", "sink"),
        (b"Name,InChI\na,SMILES=C\n", "sink"),
        (b"Name,InChI\na,InChI=1S/H2O/h1H2\nb,InChI=1S/CO2/c2-1-3\n", "source"),
    ],
)
def test_invalid_csv(payload, kind):
    with pytest.raises(InputValidationError):
        validate_compound_csv(payload, kind=kind)


def test_parameter_bounds_and_diameter_order():
    assert JobParameters().max_steps == 3
    with pytest.raises(ValueError):
        JobParameters(dmin=16, dmax=2)

