import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { resolve } from 'path'

export default defineConfig({
  plugins: [react()],
  root: resolve(__dirname, 'src/renderer'),
  base: './',
  build: {
    outDir: resolve(__dirname, 'dist/renderer'),
    emptyOutDir: true,
    rollupOptions: {
      input: resolve(__dirname, 'src/renderer/index.html')
    }
  },
  resolve: {
    alias: {
      '@': resolve(__dirname, 'src'),
      '@/components': resolve(__dirname, 'src/renderer/components'),
      '@/hooks': resolve(__dirname, 'src/renderer/hooks'),
      '@/stores': resolve(__dirname, 'src/renderer/stores'),
      '@/utils': resolve(__dirname, 'src/renderer/utils')
    }
  },
  server: {
    port: 3000,
    host: 'localhost'
  }
})