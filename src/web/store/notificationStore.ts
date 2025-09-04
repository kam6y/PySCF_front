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

export const showResourceInsufficientErrorNotification = (
  errorMessage: string,
  calculationId?: string
) => {
  // リソース不足エラーメッセージを日本語で分かりやすく変換
  let title = 'PCのリソースが不足しています';
  let message = errorMessage;
  
  // エラーメッセージの内容に応じて適切な日本語メッセージを設定
  if (errorMessage.toLowerCase().includes('cpu')) {
    if (errorMessage.includes('no active calculations')) {
      title = 'PCのCPU使用率が高すぎます';
      message = 'システムのCPU使用率が制限値を超えているため、計算を開始できません。他のプログラムを終了してCPU使用率を下げてから再試行してください。';
    } else {
      title = 'CPU使用率の制限に達しています';
      message = 'CPU使用率の制限により計算を開始できません。実行中の計算が完了するまでお待ちください。';
    }
  } else if (errorMessage.toLowerCase().includes('memory')) {
    if (errorMessage.includes('no active calculations')) {
      title = 'PCのメモリ使用量が多すぎます';
      message = 'システムのメモリ使用量が制限値を超えているため、計算を開始できません。他のプログラムを終了してメモリを確保してから再試行してください。';
    } else {
      title = 'メモリ使用量の制限に達しています';
      message = 'メモリ使用量の制限により計算を開始できません。実行中の計算が完了するまでお待ちください。';
    }
  }
  
  return useNotificationStore.getState().addNotification({
    type: 'error',
    title,
    message,
    autoClose: false, // 手動で閉じるまで表示
    duration: 0,
    calculationId,
    clickable: !!calculationId, // 計算IDがある場合はクリック可能
  });
};
