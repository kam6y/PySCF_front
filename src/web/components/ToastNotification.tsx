import React, { useEffect, useState } from 'react';
import styles from './ToastNotification.module.css';
import { Notification } from '../store/notificationStore';

interface ToastNotificationProps {
  notification: Notification;
  onClose: (id: string) => void;
  onNavigate?: (calculationId: string) => void;
}

// エラーコードに基づいて表示テキストを変換する関数
const getDisplayContent = (notification: Notification) => {
  if (!notification.errorCode) {
    return { title: notification.title, message: notification.message };
  }

  switch (notification.errorCode) {
    case 'CPU_INSUFFICIENT_SYSTEM':
      return {
        title: 'High CPU Usage',
        message: 'System CPU usage exceeds the limit. Please close other programs to reduce CPU usage before retrying.'
      };
    case 'CPU_INSUFFICIENT_LIMIT':
      return {
        title: 'CPU Limit Reached',
        message: 'Cannot start calculation due to CPU usage limit. Please wait for running calculations to complete.'
      };
    case 'MEMORY_INSUFFICIENT_SYSTEM':
      return {
        title: 'High Memory Usage',
        message: 'System memory usage exceeds the limit. Please close other programs to free up memory before retrying.'
      };
    case 'MEMORY_INSUFFICIENT_LIMIT':
      return {
        title: 'Memory Limit Reached',
        message: 'Cannot start calculation due to memory usage limit. Please wait for running calculations to complete.'
      };
    case 'RESOURCE_INSUFFICIENT':
      return {
        title: 'Insufficient Resources',
        message: notification.message || 'Insufficient system resources to start calculation.'
      };
    default:
      return { title: notification.title, message: notification.message };
  }
};

export const ToastNotification: React.FC<ToastNotificationProps> = ({
  notification,
  onClose,
  onNavigate,
}) => {
  const [isVisible, setIsVisible] = useState(false);
  const [isExiting, setIsExiting] = useState(false);
  
  const displayContent = getDisplayContent(notification);

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
        <div className={styles.toastTitle}>{displayContent.title}</div>
        {displayContent.message && (
          <div className={styles.toastMessage}>{displayContent.message}</div>
        )}
      </div>
      <button
        className={styles.toastCloseButton}
        onClick={e => {
          e.stopPropagation(); // 親のクリックイベントを阻止
          handleClose();
        }}
        aria-label="Close notification"
      >
        ×
      </button>
    </div>
  );
};
