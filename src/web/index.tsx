// src/web/index.tsx

import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { setApiBaseUrl } from "./apiClient";
import { App } from "./App";

const root = createRoot(document.getElementById("root") as Element);

// Listen for the port from the main process
window.electronAPI.onSetFlaskPort((port) => {
  console.log('Received port from main process:', port);
  setApiBaseUrl(port);
  
  // Render the app only after the port is set
  root.render(
    <StrictMode>
      <App />
    </StrictMode>
  );
});