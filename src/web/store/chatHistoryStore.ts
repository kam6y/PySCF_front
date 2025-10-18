import { create } from 'zustand';

interface ChatHistoryState {
  // 選択中のセッションID
  activeSessionId: string | null;

  // アクション
  setActiveSessionId: (sessionId: string | null) => void;
  clearActiveSession: () => void;
}

export const useChatHistoryStore = create<ChatHistoryState>((set) => ({
  // 初期状態
  activeSessionId: null,

  // アクティブセッションを設定
  setActiveSessionId: (sessionId: string | null) => set({ activeSessionId: sessionId }),

  // アクティブセッションをクリア
  clearActiveSession: () => set({ activeSessionId: null }),
}));
