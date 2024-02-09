import pytest

from pgqc.kmerlib import build_kmers

@pytest.mark.parametrize("sequence,size,expected", [
	("ATCG", 3, ["ATC", "TCG"]),
	("ATCGA", 3, ["ATC", "TCG", "CGA"]),
	("A", 3, []),
])
def test_build_kmers(sequence, size, expected):
    
	assert build_kmers(sequence, size) == expected


# def test_build_kmers_fails():

# 	with pytest.raises(ValueError):
# 		build_kmers("A", 3)