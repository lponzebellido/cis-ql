const { app, BrowserWindow, ipcMain, dialog } = require('electron');
const path = require('path');
const { spawn } = require('child_process');

let backendProcess;
let frontendProcess;
let mainWindow;
let frontendReady = false;
let detectedPort = 5173;

function createWindow(port) {
  mainWindow = new BrowserWindow({
    width: 1400,
    height: 900,
    minWidth: 900,
    minHeight: 600,
    titleBarStyle: 'hiddenInset',
    trafficLightPosition: { x: 12, y: 12 },
    backgroundColor: '#0d1117',
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      nodeIntegration: false
    }
  });

  mainWindow.loadURL(`http://localhost:${port}`);
}

app.whenReady().then(() => {
  backendProcess = spawn('node', ['backend/server.js'], {
    cwd: __dirname,
    stdio: 'pipe'
  });
  backendProcess.stdout.on('data', (data) => console.log(`[Backend] ${data}`));
  backendProcess.stderr.on('data', (data) => console.error(`[Backend] ${data}`));

  frontendProcess = spawn('npm', ['run', 'dev'], {
    cwd: path.join(__dirname, 'frontend'),
    stdio: 'pipe',
    shell: true
  });

  let stdoutBuffer = '';
  const timeoutId = setTimeout(() => {
    if (!frontendReady) {
      console.error('Frontend start timed out after 20 seconds');
      app.quit();
    }
  }, 20000);

  frontendProcess.stdout.on('data', (data) => {
    const str = data.toString();
    console.log(`[Frontend] ${str}`);

    if (!frontendReady) {
      stdoutBuffer += str;
      const match = stdoutBuffer.match(/http:\/\/localhost:(\d+)/);
      if (match) {
        frontendReady = true;
        clearTimeout(timeoutId);
        detectedPort = parseInt(match[1], 10);
        createWindow(detectedPort);
      }
    }
  });

  frontendProcess.stderr.on('data', (data) => {
    console.error(`[Frontend Error] ${data}`);
  });

  ipcMain.handle('dialog:openFolder', async () => {
    const { canceled, filePaths } = await dialog.showOpenDialog({
      title: 'Open Workspace Folder',
      properties: ['openDirectory']
    });
    if (canceled) return null;
    return filePaths[0];
  });

  app.on('activate', () => {
    if (BrowserWindow.getAllWindows().length === 0 && frontendReady) {
      createWindow(detectedPort);
    }
  });
});

app.on('window-all-closed', () => {
  if (process.platform !== 'darwin') {
    app.quit();
  }
});

app.on('quit', () => {
  if (backendProcess) backendProcess.kill();
  if (frontendProcess) frontendProcess.kill();
});
