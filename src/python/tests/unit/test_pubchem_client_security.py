"""
Regression tests for PubChem client input redaction.

Ensures that raw query strings and invalid CID values are NOT reflected
in error messages, preventing potential leakage of proprietary search
terms via logs, error strings, or crash reports.
"""

import pytest

from pubchem.client import PubChemClient, PubChemError, PubChemNotFoundError


# Sentinel values that should never appear in any error message.
PROPRIETARY_QUERY_SENTINEL = "SECRET_COMPOUND_NAME_%%%_PRIVATE"
INVALID_CID_SENTINEL = "not_a_number_SECRET_CID_%%%"


def test_invalid_cid_error_does_not_reflect_input():
    """
    GIVEN a non-numeric CID string containing a proprietary sentinel
    WHEN _find_cid is called with search_type='cid'
    THEN the raised PubChemError must NOT contain the sentinel value
    """
    client = PubChemClient()

    with pytest.raises(PubChemError) as exc_info:
        client._find_cid(INVALID_CID_SENTINEL, "cid")

    error_message = str(exc_info.value)
    assert (
        INVALID_CID_SENTINEL not in error_message
    ), f"Raw CID input was reflected in error message: {error_message}"


def test_invalid_cid_error_is_generic():
    """
    GIVEN a non-numeric CID string
    WHEN _find_cid is called with search_type='cid'
    THEN the error message should be the generic 'Invalid CID format'
    """
    client = PubChemClient()

    with pytest.raises(PubChemError, match="Invalid CID format"):
        client._find_cid("abc_not_a_number", "cid")


def test_invalid_cid_error_suppresses_chain():
    """
    GIVEN a non-numeric CID string
    WHEN _find_cid raises PubChemError
    THEN __suppress_context__ must be True (from None applied)
    """
    client = PubChemClient()

    with pytest.raises(PubChemError) as exc_info:
        client._find_cid(INVALID_CID_SENTINEL, "cid")

    assert exc_info.value.__suppress_context__ is True, (
        'PubChemError should suppress chained ValueError containing raw CID'
    )


def test_network_error_does_not_leak_query_in_error_or_log(mocker, caplog):
    """
    GIVEN a name-search query containing a proprietary sentinel
    WHEN _find_cid builds a URL embedding the sentinel and a network error occurs
    THEN neither the raised PubChemError nor the log output must contain the sentinel
    """
    import requests.exceptions
    import logging

    sentinel = PROPRIETARY_QUERY_SENTINEL
    fake_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{sentinel}/cids/JSON'
    simulated_error = requests.exceptions.ConnectionError(
        f'HTTPConnectionPool: Max retries exceeded with url: {fake_url}'
    )

    client = PubChemClient()
    mocker.patch.object(client.session, 'get', side_effect=simulated_error)

    with caplog.at_level(logging.ERROR, logger='pubchem.client'):
        with pytest.raises(PubChemError) as exc_info:
            client._find_cid(sentinel, 'name')

    error_message = str(exc_info.value)
    assert sentinel not in error_message, (
        f'Raw query leaked into PubChemError message: {error_message}'
    )

    for record in caplog.records:
        assert sentinel not in record.getMessage(), (
            f'Raw query leaked into log: {record.getMessage()}'
        )


def test_network_error_suppresses_chain(mocker):
    """
    GIVEN a network error during a name-search
    WHEN _make_request re-raises as PubChemError
    THEN __suppress_context__ must be True (from None applied)
    """
    import requests.exceptions

    sentinel = PROPRIETARY_QUERY_SENTINEL
    fake_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{sentinel}/cids/JSON'
    simulated_error = requests.exceptions.ConnectionError(
        f'HTTPConnectionPool: Max retries exceeded with url: {fake_url}'
    )

    client = PubChemClient()
    mocker.patch.object(client.session, 'get', side_effect=simulated_error)

    with pytest.raises(PubChemError) as exc_info:
        client._find_cid(sentinel, 'name')

    assert exc_info.value.__suppress_context__ is True, (
        'PubChemError should suppress chained ConnectionError containing URL'
    )


def test_not_found_error_does_not_leak_url_query(mocker):
    """
    GIVEN a URL that embeds a proprietary query sentinel
    WHEN _make_request receives a 404 response
    THEN the raised PubChemNotFoundError must NOT contain the sentinel
    """
    sentinel = PROPRIETARY_QUERY_SENTINEL
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{sentinel}/cids/JSON'

    client = PubChemClient()
    mock_response = mocker.MagicMock()
    mock_response.status_code = 404
    mocker.patch.object(client.session, 'get', return_value=mock_response)

    with pytest.raises(PubChemNotFoundError) as exc_info:
        client._make_request(url)

    error_message = str(exc_info.value)
    assert sentinel not in error_message, (
        f'Raw query leaked into PubChemNotFoundError message: {error_message}'
    )


def test_processing_error_does_not_leak_in_log_or_chain(mocker, caplog):
    """
    GIVEN a ValueError containing a sentinel during compound processing
    WHEN search_compound's except handler catches it
    THEN the sentinel must NOT appear in the raised error, the log, or the chain
    """
    import logging

    sentinel = PROPRIETARY_QUERY_SENTINEL

    client = PubChemClient()
    mocker.patch.object(client, '_find_cid', return_value=12345)
    mocker.patch.object(
        client, '_get_compound_properties',
        side_effect=ValueError(f'unexpected value: {sentinel}')
    )

    with caplog.at_level(logging.ERROR, logger='pubchem.client'):
        with pytest.raises(PubChemError) as exc_info:
            client.search_compound(sentinel, 'name')

    # Error message must be generic
    assert sentinel not in str(exc_info.value), (
        f'Sentinel in PubChemError message: {str(exc_info.value)}'
    )

    # Chain must be suppressed
    assert exc_info.value.__suppress_context__ is True, (
        'PubChemError should suppress chained ValueError containing sentinel'
    )

    # Log must not contain sentinel
    for record in caplog.records:
        assert sentinel not in record.getMessage(), (
            f'Sentinel leaked into log: {record.getMessage()}'
        )
