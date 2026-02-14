import { ApiError } from '../api/core';
import {
  showErrorNotification,
  showResourceInsufficientErrorNotification,
} from '../store/notificationStore';
import { isResourceInsufficientError } from './errorClassifier';

/**
 * Global error handler for the application.
 * Handles ApiError, Error, and unknown error types.
 *
 * @param error The error object to handle
 * @param context Optional context string to prepend to the error message (e.g. "Failed to load calculations")
 */
export const handleError = (error: unknown, context?: string) => {
  console.error('[GlobalErrorHandler]', context ? `${context}:` : '', error);

  // Default error details
  let title = 'An error occurred';
  let message = 'An unexpected error occurred. Please try again.';
  let calculationId: string | undefined = undefined;

  // Handle ApiError
  if (error instanceof ApiError) {
    calculationId = error.response?.id; // Try to extract calculation ID if available in response

    if (error.isNetworkError) {
      title = 'Network Error';
      message =
        'A network connection error occurred. Please check your internet connection.';
    } else if (error.status === 404) {
      title = 'Not Found';
      message = 'The requested resource was not found.';
    } else if (error.status === 400) {
      title = 'Invalid Request';
      message = 'There is an issue with the request. Please check your input.';
    } else if (error.status === 401) {
      title = 'Authentication Required';
      message = 'Your session has expired. Please log in again.';
    } else if (error.status === 403) {
      title = 'Access Denied';
      message = 'You do not have permission to access this resource.';
    } else if (error.status === 503) {
      title = 'Service Unavailable';
      message =
        'The server is temporarily unavailable. Please try again later.';
    } else if (error.status >= 500) {
      title = 'Server Error';
      message = 'A server error occurred. Please contact the administrator.';
    } else {
      // Use the message from the error if available, otherwise default
      message = error.message || message;
    }
  } else if (error instanceof Error) {
    // Standard Error object
    message = error.message;

    if (isResourceInsufficientError(message)) {
      showResourceInsufficientErrorNotification(message, calculationId);
      return;
    }
  } else if (typeof error === 'string') {
    message = error;
  }

  // Prepend context if provided
  if (context) {
    message = `${context}: ${message}`;
  }

  showErrorNotification(title, message, calculationId);
};
