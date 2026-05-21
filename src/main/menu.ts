import { app, Menu, shell, dialog } from 'electron';
import path from 'node:path';

/**
 * Aboutダイアログを表示する関数
 */
export const showAboutDialog = (): void => {
  // In dev: dist/main/menu.js -> ../../package.json
  const packageJson = require('../../package.json');

  const aboutMessage = `${packageJson.description}

Version: ${packageJson.version}

Key Features:
• 3D molecular structure visualization
• Quantum chemistry calculations (DFT, HF, MP2, CCSD, TDDFT)
• Geometry optimization & vibrational frequency analysis
• AI-powered multi-agent system
• PubChem/SMILES structure retrieval

Tech Stack:
Frontend: Electron, React, TypeScript
Backend: Python, FastAPI, PySCF, RDKit

License: ${packageJson.license}
Author: ${packageJson.author}

GitHub: ${packageJson.homepage}`;

  dialog.showMessageBox({
    type: 'info',
    title: 'About PySCF_front',
    message: 'PySCF_front',
    detail: aboutMessage,
    buttons: ['OK'],
    icon: app.isPackaged
      ? path.join(
          process.resourcesPath,
          'src',
          'assets',
          'icon',
          'mac',
          'Pyscf_front.icns'
        )
      : path.join(
          __dirname,
          '..',
          'src',
          'assets',
          'icon',
          'mac',
          'Pyscf_front.icns'
        ),
  });
};

/**
 * アプリケーションメニューを作成する関数
 */
export const createApplicationMenu = (): void => {
  const isMac = process.platform === 'darwin';

  const template: Electron.MenuItemConstructorOptions[] = [
    // macOSの場合、最初にアプリケーションメニューを追加
    ...(isMac
      ? [
          {
            label: app.name,
            submenu: [
              {
                label: 'About PySCF_front',
                click: showAboutDialog,
              },
              { type: 'separator' as const },
              { role: 'services' as const },
              { type: 'separator' as const },
              { role: 'hide' as const },
              { role: 'hideOthers' as const },
              { role: 'unhide' as const },
              { type: 'separator' as const },
              { role: 'quit' as const },
            ],
          },
        ]
      : []),
    // Editメニュー
    {
      label: 'Edit',
      submenu: [
        { role: 'undo' as const },
        { role: 'redo' as const },
        { type: 'separator' as const },
        { role: 'cut' as const },
        { role: 'copy' as const },
        { role: 'paste' as const },
        ...(isMac
          ? [
              { role: 'pasteAndMatchStyle' as const },
              { role: 'delete' as const },
              { role: 'selectAll' as const },
            ]
          : [
              { role: 'delete' as const },
              { type: 'separator' as const },
              { role: 'selectAll' as const },
            ]),
      ],
    },
    // Viewメニュー
    {
      label: 'View',
      submenu: [
        { role: 'reload' as const },
        { role: 'forceReload' as const },
        { role: 'toggleDevTools' as const },
        { type: 'separator' as const },
        { role: 'resetZoom' as const },
        { role: 'zoomIn' as const },
        { role: 'zoomOut' as const },
        { type: 'separator' as const },
        { role: 'togglefullscreen' as const },
      ],
    },
    // Windowメニュー
    {
      label: 'Window',
      submenu: [
        { role: 'minimize' as const },
        { role: 'zoom' as const },
        ...(isMac
          ? [
              { type: 'separator' as const },
              { role: 'front' as const },
              { type: 'separator' as const },
              { role: 'window' as const },
            ]
          : [{ role: 'close' as const }]),
      ],
    },
    // Helpメニュー
    {
      role: 'help' as const,
      submenu: [
        {
          label: 'Learn More',
          click: async () => {
            const packageJson = require('../../package.json');
            await shell.openExternal(packageJson.homepage);
          },
        },
        {
          label: 'Report Issue',
          click: async () => {
            const packageJson = require('../../package.json');
            await shell.openExternal(packageJson.bugs.url);
          },
        },
      ],
    },
  ];

  const menu = Menu.buildFromTemplate(template);
  Menu.setApplicationMenu(menu);
};
