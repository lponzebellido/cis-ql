const express = require('express');
const cors = require('cors');
const { exec } = require('child_process');
const fs = require('fs');
const path = require('path');

const app = express();
app.use(cors());
app.use(express.json());

const PORT = 3001;


const CISQL_BIN = path.resolve(__dirname, '../../cisql');
let WORKING_DIR = null;

app.get('/api/fs/workspace', (req, res) => {
  res.json({ workspace: WORKING_DIR });
});

app.post('/api/fs/workspace', (req, res) => {
  const { path: newPath } = req.body;
  if (!newPath || !fs.existsSync(newPath)) {
    return res.status(400).json({ error: 'Invalid path' });
  }
  WORKING_DIR = newPath;
  res.json({ success: true, workspace: WORKING_DIR });
});

app.post('/api/execute', (req, res) => {
  const { code, targetFile } = req.body;
  
  if (!WORKING_DIR) {
    return res.status(400).json({ error: 'No workspace opened' });
  }
  
  if (!code) {
    return res.status(400).json({ error: 'No code provided' });
  }

  
  const tempFile = path.join(WORKING_DIR, 'temp_script.cql');
  const resultsFile = path.join(WORKING_DIR, '.cisql_results.json');
  
  
  if (fs.existsSync(resultsFile)) {
    fs.unlinkSync(resultsFile);
  }

  fs.writeFileSync(tempFile, code);

  
  exec(`"${CISQL_BIN}" temp_script.cql`, { cwd: WORKING_DIR }, (error, stdout, stderr) => {
    
    let parsedResults = null;
    if (fs.existsSync(resultsFile)) {
      try {
        const rawResults = fs.readFileSync(resultsFile, 'utf8');
        parsedResults = JSON.parse(rawResults);
        fs.unlinkSync(resultsFile);
      } catch (e) {
        console.error('Error parsing results JSON:', e);
      }
    }

    if (fs.existsSync(tempFile)) {
      fs.unlinkSync(tempFile);
    }

    res.json({
      stdout: stdout,
      stderr: stderr,
      error: error ? error.message : null,
      results: parsedResults
    });
  });
});

app.get('/api/fs/list', (req, res) => {
  if (!WORKING_DIR) return res.status(400).json({ error: 'No workspace opened' });
  const dirPath = req.query.path ? path.join(WORKING_DIR, req.query.path) : WORKING_DIR;
  
  
  if (!dirPath.startsWith(WORKING_DIR)) {
    return res.status(403).json({ error: 'Forbidden' });
  }

  try {
    const files = fs.readdirSync(dirPath, { withFileTypes: true });
    
    const items = files
      .filter(f => !f.name.startsWith('.') && f.name !== 'node_modules' && f.name !== 'cisql-studio')
      .map(f => ({
        name: f.name,
        isDirectory: f.isDirectory(),
        path: path.relative(WORKING_DIR, path.join(dirPath, f.name))
      }))
      .sort((a, b) => {
        if (a.isDirectory === b.isDirectory) return a.name.localeCompare(b.name);
        return a.isDirectory ? -1 : 1;
      });
      
    res.json(items);
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

app.get('/api/fs/read', (req, res) => {
  if (!WORKING_DIR) return res.status(400).json({ error: 'No workspace opened' });
  if (!req.query.path) return res.status(400).json({ error: 'No path provided' });
  
  const targetPath = path.join(WORKING_DIR, req.query.path);
  if (!targetPath.startsWith(WORKING_DIR)) {
    return res.status(403).json({ error: 'Forbidden' });
  }

  try {
    const content = fs.readFileSync(targetPath, 'utf8');
    res.json({ content });
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

app.post('/api/fs/write', (req, res) => {
  if (!WORKING_DIR) return res.status(400).json({ error: 'No workspace opened' });
  const { path: relPath, content } = req.body;
  if (!relPath) return res.status(400).json({ error: 'No path provided' });
  
  const targetPath = path.join(WORKING_DIR, relPath);
  if (!targetPath.startsWith(WORKING_DIR)) {
    return res.status(403).json({ error: 'Forbidden' });
  }

  try {
    fs.writeFileSync(targetPath, content, 'utf8');
    res.json({ success: true });
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

app.post('/api/fs/create-file', (req, res) => {
  if (!WORKING_DIR) return res.status(400).json({ error: 'No workspace opened' });
  const { path: relPath } = req.body;
  if (!relPath) return res.status(400).json({ error: 'No path provided' });
  
  const targetPath = path.join(WORKING_DIR, relPath);
  if (!targetPath.startsWith(WORKING_DIR)) return res.status(403).json({ error: 'Forbidden' });

  try {
    if (fs.existsSync(targetPath)) {
      return res.status(400).json({ error: 'File already exists' });
    }
    fs.writeFileSync(targetPath, '', 'utf8');
    res.json({ success: true });
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

app.post('/api/fs/create-dir', (req, res) => {
  if (!WORKING_DIR) return res.status(400).json({ error: 'No workspace opened' });
  const { path: relPath } = req.body;
  if (!relPath) return res.status(400).json({ error: 'No path provided' });
  
  const targetPath = path.join(WORKING_DIR, relPath);
  if (!targetPath.startsWith(WORKING_DIR)) return res.status(403).json({ error: 'Forbidden' });

  try {
    if (fs.existsSync(targetPath)) {
      return res.status(400).json({ error: 'Directory already exists' });
    }
    fs.mkdirSync(targetPath, { recursive: true });
    res.json({ success: true });
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

app.listen(PORT, () => {
  console.log(`Cis-QL backend server running on http://localhost:${PORT}`);
});
