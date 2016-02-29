import pytest


@pytest.mark.django_db
def test_amplifiedcontent(amplifiedcontent):
    assert amplifiedcontent.name == 'human amplified content'
