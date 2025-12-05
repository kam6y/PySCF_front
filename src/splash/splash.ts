/**
 * スプラッシュウィンドウのRenderer Process用スクリプト
 * IPC通信を受信してUIを更新する
 */

// CSSのインポート
import './splash.css';

// DOM要素の取得
const stageMessageEl = document.getElementById('stageMessage') as HTMLDivElement;
const detailMessageEl = document.getElementById('detailMessage') as HTMLDivElement;
const errorContainerEl = document.getElementById('errorContainer') as HTMLDivElement;
const errorMessageEl = document.getElementById('errorMessage') as HTMLDivElement;
const splashContainerEl = document.querySelector('.splash-container') as HTMLDivElement;

/**
 * 進捗状態を更新
 */
const updateStatus = (update: { stage: string; message: string; retryCount?: number }) => {
  console.log('Splash: updateStatus', update);

  // ステージメッセージを更新
  if (stageMessageEl) {
    stageMessageEl.textContent = update.message;
  }

  // リトライカウントがある場合は詳細メッセージに表示
  if (detailMessageEl) {
    if (update.retryCount !== undefined) {
      detailMessageEl.textContent = `Attempt ${update.retryCount}`;
    } else {
      detailMessageEl.textContent = '';
    }
  }
};

/**
 * エラーを表示
 */
const showError = (message: string) => {
  console.error('Splash: showError', message);

  // エラーモードに切り替え
  if (splashContainerEl) {
    splashContainerEl.classList.add('error-mode');
  }

  // エラーメッセージを表示
  if (errorMessageEl) {
    errorMessageEl.textContent = message;
  }

  if (errorContainerEl) {
    errorContainerEl.style.display = 'flex';
  }
};

/**
 * スプラッシュをクローズ
 */
const closeSplash = () => {
  console.log('Splash: closeSplash');

  // フェードアウトアニメーションを適用
  if (splashContainerEl) {
    splashContainerEl.classList.add('fade-out');
  }
};

// IPCイベントリスナーを登録
if (window.splashAPI) {
  // 進捗状態の更新
  window.splashAPI.onUpdateStatus(updateStatus);

  // エラー表示
  window.splashAPI.onShowError(showError);

  // クローズ指示
  window.splashAPI.onClose(closeSplash);

  console.log('Splash: IPC listeners registered');
} else {
  console.error('Splash: splashAPI is not available');
}
