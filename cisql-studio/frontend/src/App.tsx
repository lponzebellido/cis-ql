import { useState, useEffect, useRef } from 'react';
import Editor, { useMonaco } from '@monaco-editor/react';
import Split from 'react-split';
import { cqlLanguageDef } from './cql-monarch';
import { TrackViewer } from './TrackViewer';
import { FileExplorer } from './FileExplorer';
import { SequenceViewer } from './SequenceViewer';
import './App.css';

const EXT_COLORS: Record<string, string> = {
  cql: '#3fb950',
  fasta: '#58a6ff',
  fa: '#58a6ff',
  gff3: '#d29922',
  gff: '#d29922',
  pwm: '#bc8cff',
  txt: '#8b949e',
  md: '#8b949e',
};

function getExtColor(filename: string): string {
  const ext = filename.split('.').pop()?.toLowerCase() || '';
  return EXT_COLORS[ext] || '#8b949e';
}

function App() {
  const [workspacePath, setWorkspacePath] = useState<string | null>(null);
  const [code, setCode] = useState('');
  const [activeFile, setActiveFile] = useState<string | null>(null);
  const [stdout, setStdout] = useState('Welcome to Cis-QL Studio.');
  const [results, setResults] = useState<Record<string, any[]>>({});
  const [gcProfiles, setGcProfiles] = useState<Record<string, any>>({});
  const [isRunning, setIsRunning] = useState(false);
  const [isSaving, setIsSaving] = useState(false);
  const [activeVisTab, setActiveVisTab] = useState<'track' | 'sequence'>('track');
  const [menuOpen, setMenuOpen] = useState<string | null>(null);
  const [isDarkMode, setIsDarkMode] = useState(true);
  const [highlightedRegion, setHighlightedRegion] = useState<any>(null);
  const monaco = useMonaco();
  const editorRef = useRef<any>(null);

  useEffect(() => {
    fetch('http://localhost:3001/api/fs/workspace')
      .then(res => res.json())
      .then(data => {
        if (data.workspace) setWorkspacePath(data.workspace);
      })
      .catch(() => setStdout('Failed to connect to backend.'));
  }, []);

  useEffect(() => {
    const mq = window.matchMedia('(prefers-color-scheme: dark)');
    setIsDarkMode(mq.matches);
    const handler = (e: MediaQueryListEvent) => setIsDarkMode(e.matches);
    mq.addEventListener('change', handler);
    return () => mq.removeEventListener('change', handler);
  }, []);

  useEffect(() => {
    if (monaco) {
      monaco.languages.register({ id: 'cql' });
      monaco.languages.setMonarchTokensProvider('cql', cqlLanguageDef as any);

      monaco.editor.defineTheme('cql-dark', {
        base: 'vs-dark',
        inherit: true,
        rules: [
          { token: 'keyword', foreground: '79c0ff', fontStyle: 'bold' },
          { token: 'type.identifier', foreground: '7ee787' },
          { token: 'string', foreground: 'a5d6ff' },
          { token: 'number', foreground: '79c0ff' },
          { token: 'comment', foreground: '8b949e' },
        ],
        colors: {
          'editor.background': '#0d1117',
          'editor.lineHighlightBackground': '#161b22',
          'editorLineNumber.foreground': '#484f58',
          'editorLineNumber.activeForeground': '#e6edf3',
        }
      });

      monaco.editor.defineTheme('cql-light', {
        base: 'vs',
        inherit: true,
        rules: [
          { token: 'keyword', foreground: '0550ae', fontStyle: 'bold' },
          { token: 'type.identifier', foreground: '116329' },
          { token: 'string', foreground: '0a3069' },
          { token: 'number', foreground: '0550ae' },
          { token: 'comment', foreground: '6e7781' },
        ],
        colors: {
          'editor.background': '#ffffff',
          'editor.lineHighlightBackground': '#f6f8fa',
          'editorLineNumber.foreground': '#8b949e',
          'editorLineNumber.activeForeground': '#1f2328',
        }
      });
    }
  }, [monaco]);

  const handleEditorDidMount = (editor: any, monacoInstance: any) => {
    editorRef.current = editor;
    editor.addCommand(monacoInstance.KeyMod.CtrlCmd | monacoInstance.KeyCode.KeyS, () => {
      handleSave();
    });
  };

  const handleSave = async () => {
    if (!activeFile) return;
    setIsSaving(true);
    try {
      const currentCode = editorRef.current.getValue();
      await fetch('http://localhost:3001/api/fs/write', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ path: activeFile, content: currentCode })
      });
      setStdout(`Saved ${activeFile}`);
    } catch (err) {
      setStdout(`Error saving file: ${err}`);
    } finally {
      setIsSaving(false);
    }
  };

  const handleExecute = async () => {
    setIsRunning(true);
    setStdout('Executing...\n');
    setResults({});
    try {
      const response = await fetch('http://localhost:3001/api/execute', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ code })
      });
      const data = await response.json();
      setStdout(data.stdout + (data.stderr ? '\n' + data.stderr : ''));
      if (data.results && data.results.resultSets) {
        setResults(data.results.resultSets);
      }
      if (data.results && data.results.gcProfiles) {
        setGcProfiles(data.results.gcProfiles);
      } else {
        setGcProfiles({});
      }
    } catch (err) {
      setStdout('Error connecting to the backend server.');
    } finally {
      setIsRunning(false);
    }
  };

  const onFileSelect = (path: string, content: string) => {
    setActiveFile(path);
    setCode(content);
  };

  const handleOpenFolder = async () => {
    let folderPath = null;
    if ((window as any).electronAPI) {
      folderPath = await (window as any).electronAPI.selectFolder();
    } else {
      folderPath = prompt("Enter the absolute path to your workspace:");
    }
    if (folderPath) {
      try {
        const res = await fetch('http://localhost:3001/api/fs/workspace', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ path: folderPath })
        });
        const data = await res.json();
        if (data.success) {
          setWorkspacePath(data.workspace);
          setActiveFile(null);
          setCode('');
        } else {
          alert(data.error);
        }
      } catch (err: any) {
        alert("Failed to set workspace: " + err.message);
      }
    }
    setMenuOpen(null);
  };

  const activeFileName = activeFile ? activeFile.split('/').pop() || '' : '';

  return (
    <div className="app-container" onClick={() => setMenuOpen(null)}>
      <div className="titlebar-drag">
        <div className="menu-brand">Cis-QL Studio</div>
        <div className="menu-items">
          <div className="menu-item-wrapper">
            <div className="menu-item" onClick={(e) => { e.stopPropagation(); setMenuOpen(menuOpen === 'file' ? null : 'file'); }}>File</div>
            {menuOpen === 'file' && (
              <div className="dropdown-menu">
                <div className="dropdown-item" onClick={handleOpenFolder}>
                  Open Folder... <span className="dropdown-shortcut">Cmd+O</span>
                </div>
                <div className="dropdown-item" onClick={() => { handleSave(); setMenuOpen(null); }}>
                  Save <span className="dropdown-shortcut">Cmd+S</span>
                </div>
              </div>
            )}
          </div>
          <div className="menu-item-wrapper">
            <div className="menu-item" onClick={(e) => { e.stopPropagation(); setMenuOpen(menuOpen === 'run' ? null : 'run'); }}>Run</div>
            {menuOpen === 'run' && (
              <div className="dropdown-menu">
                <div className="dropdown-item" onClick={() => { handleExecute(); setMenuOpen(null); }}>
                  Run Code <span className="dropdown-shortcut">Cmd+Enter</span>
                </div>
              </div>
            )}
          </div>
          <div className="menu-item-wrapper">
            <div className="menu-item" onClick={(e) => { e.stopPropagation(); setMenuOpen(menuOpen === 'terminal' ? null : 'terminal'); }}>Terminal</div>
            {menuOpen === 'terminal' && (
              <div className="dropdown-menu">
                <div className="dropdown-item" onClick={() => { setStdout(''); setMenuOpen(null); }}>
                  Clear Terminal
                </div>
              </div>
            )}
          </div>
        </div>
        <div className="top-actions">
          <span className="status-text">{isSaving ? 'Saving...' : ''}</span>
          <button
            className={`run-button ${isRunning ? 'running' : ''}`}
            onClick={handleExecute}
            disabled={isRunning || !code.trim() || !activeFile}
          >
            {isRunning ? 'Running...' : 'Run'}
          </button>
        </div>
      </div>

      <Split className="workspace" sizes={[18, 82]} minSize={140} gutterSize={3}>
        <div className="sidebar">
          <FileExplorer
            onFileSelect={onFileSelect}
            workspacePath={workspacePath}
            onSetWorkspace={handleOpenFolder}
          />
        </div>

        <Split className="main-content" direction="vertical" sizes={[72, 28]} minSize={80} gutterSize={3}>
          <Split className="editor-visualizer-split" direction="horizontal" sizes={[55, 45]} minSize={200} gutterSize={3}>
            <div className="pane editor-pane">
              <div className="pane-header">
                <div className="tab-bar">
                  {activeFile && (
                    <div className="file-tab">
                      <div className="tab-dot" style={{ background: getExtColor(activeFileName) }}></div>
                      {activeFileName}
                    </div>
                  )}
                </div>
              </div>
              <div className="pane-content">
                {activeFile ? (
                  <Editor
                    height="100%"
                    defaultLanguage="cql"
                    theme={isDarkMode ? 'cql-dark' : 'cql-light'}
                    value={code}
                    onChange={(val) => setCode(val || '')}
                    onMount={handleEditorDidMount}
                    options={{
                      minimap: { enabled: false },
                      fontSize: 13,
                      fontFamily: "'JetBrains Mono', 'Fira Code', monospace",
                      fontLigatures: true,
                      padding: { top: 16, bottom: 16 },
                      scrollBeyondLastLine: false,
                      wordWrap: 'on',
                      renderLineHighlight: 'gutter',
                      cursorBlinking: 'smooth',
                      smoothScrolling: true,
                    }}
                  />
                ) : (
                  <div className="welcome-screen">
                    <h1>Cis-QL Studio</h1>
                    <p>Cis-Regulatory Query Language</p>
                    <div className="welcome-shortcuts">
                      <div className="shortcut"><span>Open Folder</span><span>Cmd+O</span></div>
                      <div className="shortcut"><span>Save File</span><span>Cmd+S</span></div>
                      <div className="shortcut"><span>Run Code</span><span>Cmd+Enter</span></div>
                    </div>
                    {!workspacePath && (
                      <button className="welcome-open-btn" onClick={handleOpenFolder}>
                        Open Folder
                      </button>
                    )}
                  </div>
                )}
              </div>
            </div>

            <div className="pane visualizer-pane">
              <div className="pane-header">
                <div className="visualizer-tabs">
                  <div
                    className={`visualizer-tab ${activeVisTab === 'track' ? 'active' : ''}`}
                    onClick={() => setActiveVisTab('track')}
                  >
                    Track Map
                  </div>
                  <div
                    className={`visualizer-tab ${activeVisTab === 'sequence' ? 'active' : ''}`}
                    onClick={() => setActiveVisTab('sequence')}
                  >
                    Sequence
                  </div>
                </div>
              </div>
              <div className="pane-content">
                {activeVisTab === 'track' ? (
                  <TrackViewer results={results} gcProfiles={gcProfiles} onSelectRegion={(region) => {
                    setHighlightedRegion(region);
                    setActiveVisTab('sequence');
                    setTimeout(() => {
                      const el = document.getElementById('highlighted-region');
                      if (el) el.scrollIntoView({ behavior: 'smooth', block: 'center' });
                    }, 100);
                  }} />
                ) : (
                  <SequenceViewer results={results} highlightedRegion={highlightedRegion} />
                )}
              </div>
            </div>
          </Split>

          <div className="pane terminal-pane">
            <div className="pane-header">
              <div className="pane-label">Terminal</div>
            </div>
            <pre className="terminal-output">{stdout}</pre>
          </div>
        </Split>
      </Split>

      <div className="status-bar">
        <div className="status-item">{workspacePath ? workspacePath.split('/').pop() : 'No workspace'}</div>
        {activeFile && (
          <div className="status-item">{activeFile}</div>
        )}
        <div className="spacer" />
        <div className="status-item">UTF-8</div>
        <div className="status-item">Cis-QL</div>
      </div>
    </div>
  );
}

export default App;
