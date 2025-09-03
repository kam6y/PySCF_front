import React, { useEffect, useState } from 'react';
import styles from './ToastNotification.module.css';
import { Notification } from '../store/notificationStore';

interface ToastNotificationProps {
  notification: Notification;
  onClose: (id: string) => void;
  onNavigate?: (calculationId: string) => void;
}

export const ToastNotification: React.FC<ToastNotificationProps> = ({
  notification,
  onClose,
  onNavigate,
}) => {
  const [isVisible, setIsVisible] = useState(false);
  const [isExiting, setIsExiting] = useState(false);

  useEffect(() => {
    // Trigger entrance animation
    const timer = setTimeout(() => setIsVisible(true), 10);
    return () => clearTimeout(timer);
  }, []);

  const handleClose = () => {
    setIsExiting(true);
    // Wait for exit animation to complete
    setTimeout(() => {
      onClose(notification.id);
    }, 300);
  };

  const handleClick = () => {
    if (notification.clickable && notification.calculationId && onNavigate) {
      onNavigate(notification.calculationId);
    }
  };

  const getIcon = () => {
    switch (notification.type) {
      case 'error':
        return '❌';
      case 'success':
        return '✅';
      case 'info':
        return 'ℹ️';
      default:
        return 'ℹ️';
    }
  };

  const toastClassName = [
    styles.toast,
    styles[
      `toast${notification.type.charAt(0).toUpperCase()}${notification.type.slice(1)}`
    ],
    isVisible ? styles.toastVisible : '',
    isExiting ? styles.toastExiting : '',
    notification.clickable ? styles.toastClickable : '',
  ]
    .filter(Boolean)
    .join(' ');

  return (
    <div
      className={toastClassName}
      onClick={handleClick}
      style={{ cursor: notification.clickable ? 'pointer' : 'default' }}
    >
      <div className={styles.toastIcon}>{getIcon()}</div>
      <div className={styles.toastContent}>
        <div className={styles.toastTitle}>{notification.title}</div>
        {notification.message && (
          <div className={styles.toastMessage}>{notification.message}</div>
        )}
      </div>
      <button
        className={styles.toastCloseButton}
        onClick={e => {
          e.stopPropagation(); // 親のクリックイベントを阻止
          handleClose();
        }}
        aria-label="通知を閉じる"
      >
        ×
      </button>
    </div>
  );
};
