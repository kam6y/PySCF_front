import { create } from 'zustand';

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
  // 会話履歴（メモリキャッシュのみ、永続化はDBで行う）
  history: ChatHistory[];

  // エージェントステータス
  currentAgentStatus: AgentStatus;

  // アクション
  addMessage: (message: ChatHistory) => void;
  addMessages: (messages: ChatHistory[]) => void;
  updateMessage: (index: number, update: Partial<ChatHistory>) => void;
  setHistory: (history: ChatHistory[]) => void;
  clearHistory: () => void;
  setAgentStatus: (status: AgentStatus) => void;
}

// 永続化を削除し、メモリキャッシュとしてのみ使用
// 会話履歴の永続化はSQLiteデータベースで一元管理
// セッションIDの管理はchatHistoryStoreで一元化
export const useAgentStore = create<AgentState>((set, get) => ({
  // 初期状態
  history: [],
  currentAgentStatus: { agent: null, status: 'idle' },

  // メッセージを履歴に追加
  addMessage: (message: ChatHistory) =>
    set(state => ({
      history: [...state.history, message],
    })),

  // 複数のメッセージを一度に履歴に追加（状態更新の原子性を保証）
  addMessages: (messages: ChatHistory[]) =>
    set(state => ({
      history: [...state.history, ...messages],
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
}));
