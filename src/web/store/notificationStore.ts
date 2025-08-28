import { create } from 'zustand';

export interface Notification {
  id: string;
  type: 'error' | 'success' | 'info';
  title: string;
  message?: string;
  autoClose: boolean;
  duration: number;
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
export const showErrorNotification = (title: string, message?: string) => {
  return useNotificationStore.getState().addNotification({
    type: 'error',
    title,
    message,
    autoClose: false, // 無限表示（手動で閉じるまで表示）
    duration: 0, // autoCloseがfalseの場合、durationは使用されない
  });
};

export const showSuccessNotification = (title: string, message?: string) => {
  return useNotificationStore.getState().addNotification({
    type: 'success',
    title,
    message,
    autoClose: false,
    duration: 0,
  });
};

export const showInfoNotification = (title: string, message?: string) => {
  return useNotificationStore.getState().addNotification({
    type: 'info',
    title,
    message,
    autoClose: true,
    duration: 4000,
  });
};
