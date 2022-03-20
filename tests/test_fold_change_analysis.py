import pytest
import fold_change_analysis

def test_find_unique_genes():
    fruits = ["apple", "orange", "tomato"]
    vegetables = ["tomato", "celery", "spinach"]

    assert fold_change_analysis.find_unique_genes(fruits, vegetables) == ["apple", "orange", "celery", "spinach"]
