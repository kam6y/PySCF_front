const initialPort = window.electronAPI?.flaskPort;

if (!initialPort) {
  console.error('[API Client] Flask port not available. API calls will fail explicitly.');
}

// フォールバック 5000 を廃止。ポート未取得時は空文字列にして fetch が明確に失敗するようにする
let API_BASE_URL = initialPort ? `http://127.0.0.1:${initialPort}` : '';

console.log(`[API Client] Initialized with port: ${initialPort ?? 'UNAVAILABLE'}`);

export const setApiBaseUrl = (port: number) => {
  API_BASE_URL = `http://127.0.0.1:${port}`;
  console.log(`[API Client] API base URL updated to: ${API_BASE_URL}`);
};

export const getApiBaseUrl = () => API_BASE_URL;

export type ApiResponse<T> = {
  success: boolean;
  data: T;
  error?: string;
};

export class ApiError extends Error {
  public readonly status: number;
  public readonly statusText: string;
  public readonly url: string;
  public readonly response?: unknown;
  public readonly isNetworkError: boolean;

  constructor(
    message: string,
    status: number,
    statusText: string,
    url: string,
    response?: unknown,
    isNetworkError = false
  ) {
    super(message);
    this.name = 'ApiError';
    this.status = status;
    this.statusText = statusText;
    this.url = url;
    this.response = response;
    this.isNetworkError = isNetworkError;

    if (Error.captureStackTrace) {
      Error.captureStackTrace(this, ApiError);
    }
  }
}

export const request = async <T>(
  endpoint: string,
  options: RequestInit = {}
): Promise<T> => {
  const url = `${API_BASE_URL}${endpoint}`;

  try {
    const authToken = await window.electronAPI?.getAuthToken();
    const headers: HeadersInit = {
      'Content-Type': 'application/json',
      ...options.headers,
    };

    if (authToken) {
      (headers as Record<string, string>)['X-Auth-Token'] = authToken;
    }

    const response = await fetch(url, {
      headers,
      ...options,
    });

    let apiResponse: ApiResponse<T>;

    try {
      apiResponse = await response.json();
    } catch (jsonError) {
      throw new ApiError(
        'An invalid response was returned from the server.',
        response.status,
        response.statusText,
        url,
        null,
        false
      );
    }

    if (!response.ok || !apiResponse.success) {
      const errorMessage =
        apiResponse.error ||
        `HTTPエラー: ${response.status} ${response.statusText}`;

      throw new ApiError(
        errorMessage,
        response.status,
        response.statusText,
        url,
        apiResponse,
        false
      );
    }

    return apiResponse.data;
  } catch (error) {
    if (error instanceof ApiError) {
      throw error;
    }

    throw new ApiError(
      'A network error occurred. Unable to connect to the server.',
      0,
      'Network Error',
      url,
      null,
      true
    );
  }
};
