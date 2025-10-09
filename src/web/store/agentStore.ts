import { create } from 'zustand';
import { persist } from 'zustand/middleware';

// 会話履歴の型定義
export type ChatHistory = {
  role: 'user' | 'model';
  parts: { text: string }[];
  isStreaming?: boolean;
  tempId?: number;
};

// エージェントステータスの型定義
export type AgentStatus = {
  agent: 'supervisor' | 'quantum_calculation_worker' | 'research_expert' | null;
  status: 'idle' | 'running' | 'responding';
};

interface AgentState {
  // 会話履歴
  history: ChatHistory[];

  // エージェントステータス
  currentAgentStatus: AgentStatus;

  // アクション
  addMessage: (message: ChatHistory) => void;
  updateMessage: (index: number, update: Partial<ChatHistory>) => void;
  setHistory: (history: ChatHistory[]) => void;
  clearHistory: () => void;
  setAgentStatus: (status: AgentStatus) => void;
}

export const useAgentStore = create<AgentState>()(
  persist(
    (set, get) => ({
      // 初期状態
      history: [],
      currentAgentStatus: { agent: null, status: 'idle' },

      // メッセージを履歴に追加
      addMessage: (message: ChatHistory) =>
        set(state => ({
          history: [...state.history, message],
        })),

      // 特定のインデックスのメッセージを更新（ストリーミング用）
      updateMessage: (index: number, update: Partial<ChatHistory>) =>
        set(state => {
          const newHistory = [...state.history];
          if (index >= 0 && index < newHistory.length) {
            newHistory[index] = { ...newHistory[index], ...update };
          }
          return { history: newHistory };
        }),

      // 履歴全体を設定（置き換え）
      setHistory: (history: ChatHistory[]) => set({ history }),

      // 会話履歴をクリア
      clearHistory: () => set({ history: [] }),

      // エージェントステータスを設定
      setAgentStatus: (status: AgentStatus) => set({ currentAgentStatus: status }),
    }),
    {
      name: 'pyscf-agent-storage', // localStorageのキー名
      // isStreaming、tempId、currentAgentStatusは永続化しない
      partialize: state => ({
        history: state.history.map(({ isStreaming, tempId, ...rest }) => rest),
      }),
    }
  )
);
