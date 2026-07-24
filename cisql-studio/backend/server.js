const express = require('express');
const cors = require('cors');
const { execFile } = require('child_process');
const fs = require('fs');
const path = require('path');

const app = express();
app.use(cors());
app.use(express.json());

const PORT = 3001;


const CISQL_BIN = path.resolve(__dirname, '../../cisql');
let WORKING_DIR = null;
let executionInProgress = false;

function isInside(root, candidate) {
  const relative = path.relative(root, candidate);
  return relative === '' ||
    (!relative.startsWith(`..${path.sep}`) && relative !== '..' &&
     !path.isAbsolute(relative));
}

function workspacePath(relativePath = '.') {
  if (!WORKING_DIR) return null;
  const root = fs.realpathSync(WORKING_DIR);
  const candidate = path.resolve(root, relativePath);
  if (!isInside(root, candidate)) return null;

  let existingAncestor = candidate;
  while (!fs.existsSync(existingAncestor) && existingAncestor !== root) {
    existingAncestor = path.dirname(existingAncestor);
  }
  if (!isInside(root, fs.realpathSync(existingAncestor))) return null;
  if (fs.existsSync(candidate) &&
      !isInside(root, fs.realpathSync(candidate))) return null;
  return candidate;
}

app.get('/api/fs/workspace', (req, res) => {
  res.json({ workspace: WORKING_DIR });
});

app.post('/api/fs/workspace', (req, res) => {
  const { path: newPath } = req.body;
  if (!newPath || !fs.existsSync(newPath) ||
      !fs.statSync(newPath).isDirectory()) {
    return res.status(400).json({ error: 'Invalid path' });
  }
  WORKING_DIR = fs.realpathSync(newPath);
  res.json({ success: true, workspace: WORKING_DIR });
});

app.post('/api/execute', (req, res) => {
  const { code } = req.body;
  
  if (!WORKING_DIR) {
    return res.status(400).json({ error: 'No workspace opened' });
  }
  
  if (!code) {
    return res.status(400).json({ error: 'No code provided' });
  }
  if (executionInProgress) {
    return res.status(409).json({ error: 'A query is already running' });
  }

  executionInProgress = true;
  const tempFile = path.join(
    WORKING_DIR, `.cisql-studio-${process.pid}-${Date.now()}.cql`
  );
  const resultsFile = path.join(WORKING_DIR, '.cisql_results.json');
  
  
  try {
    if (fs.existsSync(resultsFile)) {
      fs.unlinkSync(resultsFile);
    }
    fs.writeFileSync(tempFile, code);
  } catch (error) {
    executionInProgress = false;
    return res.status(500).json({ error: error.message });
  }

  
  execFile(CISQL_BIN, [path.basename(tempFile)], { cwd: WORKING_DIR },
    (error, stdout, stderr) => {
    
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
    executionInProgress = false;

    res.json({
      stdout: stdout,
      stderr: stderr,
      error: error ? error.message : null,
      exitCode: error && typeof error.code === 'number' ? error.code : 0,
      results: parsedResults
    });
  });
});

app.get('/api/fs/list', (req, res) => {
  if (!WORKING_DIR) return res.status(400).json({ error: 'No workspace opened' });
  const dirPath = workspacePath(req.query.path || '.');
  if (!dirPath) {
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
  
  const targetPath = workspacePath(req.query.path);
  if (!targetPath) {
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
  
  const targetPath = workspacePath(relPath);
  if (!targetPath) {
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
  
  const targetPath = workspacePath(relPath);
  if (!targetPath) return res.status(403).json({ error: 'Forbidden' });

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
  
  const targetPath = workspacePath(relPath);
  if (!targetPath) return res.status(403).json({ error: 'Forbidden' });

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

app.post('/api/fs/rename', (req, res) => {
  if (!WORKING_DIR) return res.status(400).json({ error: 'No workspace opened' });
  const { oldPath, newPath } = req.body;
  if (!oldPath || !newPath) return res.status(400).json({ error: 'Paths not provided' });

  const targetOldPath = workspacePath(oldPath);
  const targetNewPath = workspacePath(newPath);

  if (!targetOldPath || !targetNewPath ||
      targetOldPath === WORKING_DIR || targetNewPath === WORKING_DIR) {
    return res.status(403).json({ error: 'Forbidden' });
  }

  try {
    if (!fs.existsSync(targetOldPath)) {
      return res.status(404).json({ error: 'Source file does not exist' });
    }
    if (fs.existsSync(targetNewPath)) {
      return res.status(400).json({ error: 'Destination already exists' });
    }
    fs.renameSync(targetOldPath, targetNewPath);
    res.json({ success: true });
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

app.post('/api/fs/delete', (req, res) => {
  if (!WORKING_DIR) return res.status(400).json({ error: 'No workspace opened' });
  const { path: relPath } = req.body;
  if (!relPath) return res.status(400).json({ error: 'No path provided' });

  const targetPath = workspacePath(relPath);
  if (!targetPath || targetPath === WORKING_DIR) {
    return res.status(403).json({ error: 'Forbidden' });
  }

  try {
    if (!fs.existsSync(targetPath)) {
      return res.status(404).json({ error: 'File does not exist' });
    }
    const stat = fs.statSync(targetPath);
    if (stat.isDirectory()) {
      fs.rmSync(targetPath, { recursive: true, force: true });
    } else {
      fs.unlinkSync(targetPath);
    }
    res.json({ success: true });
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

app.listen(PORT, () => {
  console.log(`Cis-QL backend server running on http://localhost:${PORT}`);
});
