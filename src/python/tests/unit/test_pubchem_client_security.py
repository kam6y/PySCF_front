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


# ============================================================================
# URL Path Encoding Tests (L-001)
# ============================================================================

class TestPubChemQueryUrlEncoding:
    """Verify that user-supplied query strings are percent-encoded in URL paths.

    PubChem queries are embedded in REST URL path segments.  Characters such as
    ``/``, ``?``, ``#``, and non-ASCII must be percent-encoded so they cannot
    alter path or query semantics of the upstream request.
    """

    def test_slash_in_query_is_encoded(self, mocker):
        """
        GIVEN a query containing a forward slash
        WHEN _find_cid builds the URL
        THEN the slash must be percent-encoded (%2F) in the path
        """
        client = PubChemClient()
        captured_url = {}

        def _capture_get(url, **kwargs):
            captured_url['url'] = url
            mock_resp = mocker.MagicMock()
            mock_resp.status_code = 200
            mock_resp.json.return_value = {'IdentifierList': {'CID': [12345]}}
            return mock_resp

        mocker.patch.object(client.session, 'get', side_effect=_capture_get)
        client._find_cid('a/b', 'name')

        assert '%2F' in captured_url['url'], (
            f'Slash was not percent-encoded in URL: {captured_url["url"]}'
        )
        assert '/a/b/' not in captured_url['url'], (
            'Raw slash appeared as a path separator in URL'
        )

    def test_question_mark_in_query_is_encoded(self, mocker):
        """
        GIVEN a query containing a question mark
        WHEN _find_cid builds the URL
        THEN the question mark must be percent-encoded (%3F)
        """
        client = PubChemClient()
        captured_url = {}

        def _capture_get(url, **kwargs):
            captured_url['url'] = url
            mock_resp = mocker.MagicMock()
            mock_resp.status_code = 200
            mock_resp.json.return_value = {'IdentifierList': {'CID': [12345]}}
            return mock_resp

        mocker.patch.object(client.session, 'get', side_effect=_capture_get)
        client._find_cid('water?extra=1', 'name')

        assert '%3F' in captured_url['url'], (
            f'Question mark was not percent-encoded in URL: {captured_url["url"]}'
        )

    def test_hash_in_query_is_encoded(self, mocker):
        """
        GIVEN a query containing a hash
        WHEN _find_cid builds the URL
        THEN the hash must be percent-encoded (%23)
        """
        client = PubChemClient()
        captured_url = {}

        def _capture_get(url, **kwargs):
            captured_url['url'] = url
            mock_resp = mocker.MagicMock()
            mock_resp.status_code = 200
            mock_resp.json.return_value = {'IdentifierList': {'CID': [12345]}}
            return mock_resp

        mocker.patch.object(client.session, 'get', side_effect=_capture_get)
        client._find_cid('ethanol#fragment', 'name')

        assert '%23' in captured_url['url'], (
            f'Hash was not percent-encoded in URL: {captured_url["url"]}'
        )

    def test_space_in_query_is_encoded(self, mocker):
        """
        GIVEN a query containing a space
        WHEN _find_cid builds the URL
        THEN the space must be percent-encoded (%20)
        """
        client = PubChemClient()
        captured_url = {}

        def _capture_get(url, **kwargs):
            captured_url['url'] = url
            mock_resp = mocker.MagicMock()
            mock_resp.status_code = 200
            mock_resp.json.return_value = {'IdentifierList': {'CID': [12345]}}
            return mock_resp

        mocker.patch.object(client.session, 'get', side_effect=_capture_get)
        client._find_cid('acetic acid', 'name')

        assert '%20' in captured_url['url'], (
            f'Space was not percent-encoded in URL: {captured_url["url"]}'
        )
        assert ' ' not in captured_url['url'].split('?')[0], (
            'Raw space appeared in URL path'
        )

    def test_non_ascii_in_query_is_encoded(self, mocker):
        """
        GIVEN a query containing non-ASCII characters (e.g. CJK compound name)
        WHEN _find_cid builds the URL
        THEN the non-ASCII bytes must be percent-encoded
        """
        client = PubChemClient()
        captured_url = {}

        def _capture_get(url, **kwargs):
            captured_url['url'] = url
            mock_resp = mocker.MagicMock()
            mock_resp.status_code = 200
            mock_resp.json.return_value = {'IdentifierList': {'CID': [12345]}}
            return mock_resp

        mocker.patch.object(client.session, 'get', side_effect=_capture_get)
        client._find_cid('アスピリン', 'name')  # katakana 'aspirin'

        # Non-ASCII chars should be percent-encoded; raw katakana should not appear
        assert 'ア' not in captured_url['url'], (
            'Non-ASCII character appeared un-encoded in URL'
        )

    def test_cid_search_type_bypasses_url_encoding(self):
        """
        GIVEN search_type='cid' with a valid numeric CID
        WHEN _find_cid is called
        THEN it returns the integer directly without building a URL
        """
        client = PubChemClient()
        result = client._find_cid('12345', 'cid')
        assert result == 12345

    def test_combined_special_characters(self, mocker):
        """
        GIVEN a query combining /, ?, #, and spaces
        WHEN _find_cid builds the URL
        THEN all special characters are percent-encoded
        """
        client = PubChemClient()
        captured_url = {}

        def _capture_get(url, **kwargs):
            captured_url['url'] = url
            mock_resp = mocker.MagicMock()
            mock_resp.status_code = 200
            mock_resp.json.return_value = {'IdentifierList': {'CID': [12345]}}
            return mock_resp

        mocker.patch.object(client.session, 'get', side_effect=_capture_get)
        client._find_cid('a/b?c#d e', 'name')

        url = captured_url['url']
        assert '%2F' in url, 'Slash not encoded'
        assert '%3F' in url, 'Question mark not encoded'
        assert '%23' in url, 'Hash not encoded'
        assert '%20' in url, 'Space not encoded'


# ============================================================================
# Search Type Validation Tests (G-001)
# ============================================================================

class TestPubChemSearchTypeValidation:
    """Verify that _find_cid rejects invalid search_type values at the client boundary.

    The service layer already validates search_type, but the client method is a
    public surface and defense-in-depth requires validation at every trust boundary.
    """

    def test_invalid_search_type_raises_error(self):
        """
        GIVEN an invalid search_type string
        WHEN _find_cid is called
        THEN PubChemError is raised with a generic message
        """
        client = PubChemClient()
        with pytest.raises(PubChemError, match='Invalid search type'):
            client._find_cid('water', 'invalid_type')

    def test_path_traversal_search_type_rejected(self):
        """
        GIVEN a search_type containing path traversal characters
        WHEN _find_cid is called
        THEN PubChemError is raised (not interpolated into URL)
        """
        client = PubChemClient()
        with pytest.raises(PubChemError, match='Invalid search type'):
            client._find_cid('water', '../../../etc/passwd')

    def test_empty_search_type_rejected(self):
        """
        GIVEN an empty string search_type
        WHEN _find_cid is called
        THEN PubChemError is raised
        """
        client = PubChemClient()
        with pytest.raises(PubChemError, match='Invalid search type'):
            client._find_cid('water', '')

    @pytest.mark.parametrize('valid_type', ['name', 'cid', 'formula'])
    def test_valid_search_types_accepted(self, valid_type, mocker):
        """
        GIVEN a valid search_type
        WHEN _find_cid is called
        THEN no PubChemError about invalid search type is raised
        """
        client = PubChemClient()
        if valid_type == 'cid':
            # CID path returns int directly, no network call needed
            result = client._find_cid('12345', 'cid')
            assert result == 12345
        else:
            # Mock network call for name/formula
            mock_resp = mocker.MagicMock()
            mock_resp.status_code = 200
            mock_resp.json.return_value = {'IdentifierList': {'CID': [99999]}}
            mocker.patch.object(client.session, 'get', return_value=mock_resp)
            result = client._find_cid('test', valid_type)
            assert result == 99999
