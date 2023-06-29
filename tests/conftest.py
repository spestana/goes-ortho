import pytest
from goes_ortho import get_data

@pytest.fixture(scope="session")
def setup_session(request):
    """Setup test state for this test session."""
    print("\nDoing setup")
    get_data.download_example_data()
    def fin():
        """Teardown test state for this test session."""
        print("\nDoing teardown")
        get_data.remove_example_data()
    request.addfinalizer(fin)