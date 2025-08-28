import React from 'react';
import { useNotificationStore } from '../store/notificationStore';
import { ToastNotification } from './ToastNotification';
import styles from './ToastContainer.module.css';

export const ToastContainer: React.FC = () => {
  const { notifications, removeNotification } = useNotificationStore();

  if (notifications.length === 0) {
    return null;
  }

  return (
    <div className={styles.toastContainer}>
      {notifications.map(notification => (
        <ToastNotification
          key={notification.id}
          notification={notification}
          onClose={removeNotification}
        />
      ))}
    </div>
  );
};
