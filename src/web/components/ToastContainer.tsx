import React from 'react';
import { useNotificationStore } from '../store/notificationStore';
import { ToastNotification } from './ToastNotification';
import { useAppState } from '../hooks/useAppState';
import styles from './ToastContainer.module.css';

export const ToastContainer: React.FC = () => {
  const { notifications, removeNotification } = useNotificationStore();
  const appState = useAppState();

  const handleNavigate = (calculationId: string) => {
    // 該当する計算を選択
    appState.calculation.selectCalculation(calculationId);

    // 適切なページに遷移
    // 成功通知の場合は結果ページ、エラー通知の場合は設定ページに遷移
    const notification = notifications.find(
      n => n.calculationId === calculationId
    );
    if (notification) {
      if (notification.type === 'success') {
        appState.ui.setCurrentPage('calculation-results');
      } else if (notification.type === 'error') {
        appState.ui.setCurrentPage('calculation-settings');
      }
    }
  };

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
          onNavigate={handleNavigate}
        />
      ))}
    </div>
  );
};
