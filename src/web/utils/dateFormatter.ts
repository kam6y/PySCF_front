/**
 * Date and time formatting utilities with timezone support.
 */

/**
 * Timezone display labels with both abbreviation and full name.
 */
export const TIMEZONE_LABELS: Record<string, string> = {
  UTC: 'UTC (Coordinated Universal Time)',
  'Asia/Tokyo': 'JST (Japan Standard Time)',
  'Asia/Shanghai': 'CST (China Standard Time)',
  'Asia/Seoul': 'KST (Korea Standard Time)',
  'Asia/Singapore': 'SGT (Singapore Time)',
  'Asia/Kolkata': 'IST (India Standard Time)',
  'Australia/Sydney': 'AEST (Australian Eastern Standard Time)',
  'Europe/London': 'GMT (Greenwich Mean Time)',
  'Europe/Paris': 'CET (Central European Time)',
  'Europe/Berlin': 'CET (Central European Time)',
  'America/New_York': 'EST (Eastern Standard Time)',
  'America/Chicago': 'CST (Central Standard Time)',
  'America/Denver': 'MST (Mountain Standard Time)',
  'America/Los_Angeles': 'PST (Pacific Standard Time)',
};

/**
 * Timezone groups for UI organization.
 */
export const TIMEZONE_GROUPS = {
  standard: ['UTC'],
  asiaPacific: [
    'Asia/Tokyo',
    'Asia/Shanghai',
    'Asia/Seoul',
    'Asia/Singapore',
    'Asia/Kolkata',
    'Australia/Sydney',
  ],
  europe: ['Europe/London', 'Europe/Paris', 'Europe/Berlin'],
  americas: [
    'America/New_York',
    'America/Chicago',
    'America/Denver',
    'America/Los_Angeles',
  ],
};

/**
 * Get the display label for a timezone.
 *
 * @param timezone - IANA timezone identifier
 * @returns Display label with abbreviation and full name
 */
export function getTimezoneLabel(timezone: string): string {
  return TIMEZONE_LABELS[timezone] || timezone;
}

/**
 * Format a date-time string according to the specified timezone.
 *
 * @param dateString - ISO 8601 date-time string
 * @param timezone - IANA timezone identifier (e.g., 'UTC', 'Asia/Tokyo')
 * @returns Formatted date-time string
 *
 * @example
 * formatDateTime('2024-01-15T10:30:00Z', 'UTC')
 * // Returns: "01/15/2024, 10:30"
 *
 * formatDateTime('2024-01-15T10:30:00Z', 'Asia/Tokyo')
 * // Returns: "2024/01/15 19:30"
 */
export function formatDateTime(dateString: string, timezone: string): string {
  try {
    const date = new Date(dateString);

    // Check if date is valid
    if (isNaN(date.getTime())) {
      console.warn(`Invalid date string: ${dateString}`);
      return dateString;
    }

    // Determine locale based on timezone
    const locale = timezone === 'Asia/Tokyo' ? 'ja-JP' : 'en-US';

    // Format options
    const options: Intl.DateTimeFormatOptions = {
      timeZone: timezone,
      year: 'numeric',
      month: '2-digit',
      day: '2-digit',
      hour: '2-digit',
      minute: '2-digit',
    };

    return date.toLocaleString(locale, options);
  } catch (error) {
    console.error(`Error formatting date: ${error}`);
    return dateString;
  }
}

/**
 * Format a date-time string with seconds.
 *
 * @param dateString - ISO 8601 date-time string
 * @param timezone - IANA timezone identifier
 * @returns Formatted date-time string with seconds
 */
export function formatDateTimeWithSeconds(
  dateString: string,
  timezone: string
): string {
  try {
    const date = new Date(dateString);

    if (isNaN(date.getTime())) {
      console.warn(`Invalid date string: ${dateString}`);
      return dateString;
    }

    const locale = timezone === 'Asia/Tokyo' ? 'ja-JP' : 'en-US';

    const options: Intl.DateTimeFormatOptions = {
      timeZone: timezone,
      year: 'numeric',
      month: '2-digit',
      day: '2-digit',
      hour: '2-digit',
      minute: '2-digit',
      second: '2-digit',
    };

    return date.toLocaleString(locale, options);
  } catch (error) {
    console.error(`Error formatting date with seconds: ${error}`);
    return dateString;
  }
}

/**
 * Get the current date-time in the specified timezone.
 *
 * @param timezone - IANA timezone identifier
 * @returns Formatted current date-time string
 */
export function getCurrentDateTime(timezone: string): string {
  return formatDateTime(new Date().toISOString(), timezone);
}

/**
 * Check if a timezone is valid.
 *
 * @param timezone - IANA timezone identifier
 * @returns True if timezone is valid
 */
export function isValidTimezone(timezone: string): boolean {
  return timezone in TIMEZONE_LABELS;
}
