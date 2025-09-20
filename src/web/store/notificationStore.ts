import { create } from 'zustand';

export interface Notification {
  id: string;
  type: 'error' | 'success' | 'info';
  title: string;
  message?: string;
  autoClose: boolean;
  duration: number;
  calculationId?: string;
  clickable?: boolean;
  errorCode?: string;
}

interface NotificationState {
  notifications: Notification[];
  addNotification: (notification: Omit<Notification, 'id'>) => string;
  removeNotification: (id: string) => void;
  clearNotifications: () => void;
}

export const useNotificationStore = create<NotificationState>((set, get) => ({
  notifications: [],

  addNotification: notificationData => {
    const id = `notification-${Date.now()}-${Math.random().toString(36).substring(2, 11)}`;
    const notification: Notification = {
      id,
      ...notificationData,
    };

    set(state => ({
      notifications: [...state.notifications, notification],
    }));

    // Auto-close if specified
    if (notification.autoClose) {
      setTimeout(() => {
        get().removeNotification(id);
      }, notification.duration);
    }

    return id;
  },

  removeNotification: id => {
    set(state => ({
      notifications: state.notifications.filter(n => n.id !== id),
    }));
  },

  clearNotifications: () => {
    set({ notifications: [] });
  },
}));

// Convenience functions for common notification types
export const showErrorNotification = (
  title: string,
  message?: string,
  calculationId?: string
) => {
  return useNotificationStore.getState().addNotification({
    type: 'error',
    title,
    message,
    autoClose: false, // 無限表示（手動で閉じるまで表示）
    duration: 0, // autoCloseがfalseの場合、durationは使用されない
    calculationId,
    clickable: !!calculationId, // 計算IDがある場合はクリック可能
  });
};

export const showSuccessNotification = (
  title: string,
  message?: string,
  calculationId?: string
) => {
  return useNotificationStore.getState().addNotification({
    type: 'success',
    title,
    message,
    autoClose: false,
    duration: 0,
    calculationId,
    clickable: !!calculationId, // 計算IDがある場合はクリック可能
  });
};

export const showInfoNotification = (
  title: string,
  message?: string,
  calculationId?: string
) => {
  return useNotificationStore.getState().addNotification({
    type: 'info',
    title,
    message,
    autoClose: true,
    duration: 4000,
    calculationId,
    clickable: !!calculationId, // 計算IDがある場合はクリック可能
  });
};

// エラーコード付きエラー通知の汎用関数
export const showErrorWithCodeNotification = (
  errorCode: string,
  rawErrorMessage: string,
  title?: string,
  calculationId?: string
) => {
  return useNotificationStore.getState().addNotification({
    type: 'error',
    title: title || 'Error',
    message: rawErrorMessage,
    autoClose: false,
    duration: 0,
    calculationId,
    clickable: !!calculationId,
    errorCode,
  });
};

// 後方互換性のためのリソース不足エラー通知関数
export const showResourceInsufficientErrorNotification = (
  errorMessage: string,
  calculationId?: string
) => {
  // エラーメッセージの内容からエラーコードを判定
  let errorCode = 'RESOURCE_INSUFFICIENT';
  
  if (errorMessage.toLowerCase().includes('cpu')) {
    errorCode = errorMessage.includes('no active calculations') 
      ? 'CPU_INSUFFICIENT_SYSTEM' 
      : 'CPU_INSUFFICIENT_LIMIT';
  } else if (errorMessage.toLowerCase().includes('memory')) {
    errorCode = errorMessage.includes('no active calculations') 
      ? 'MEMORY_INSUFFICIENT_SYSTEM' 
      : 'MEMORY_INSUFFICIENT_LIMIT';
  }

  return showErrorWithCodeNotification(errorCode, errorMessage, 'Resource Insufficient', calculationId);
};
