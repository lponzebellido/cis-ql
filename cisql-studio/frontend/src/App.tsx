import { useState, useEffect, useRef } from 'react';
import Editor, { useMonaco } from '@monaco-editor/react';
import Split from 'react-split';
import { cqlLanguageDef } from './cql-monarch';
import { TrackViewer } from './TrackViewer';
import { FileExplorer } from './FileExplorer';
import { SequenceViewer } from './SequenceViewer';
import './App.css';

function App() {
  const [workspacePath, setWorkspacePath] = useState<string | null>(null);
  const [code, setCode] = useState('');
  const [activeFile, setActiveFile] = useState<string | null>(null);
  const [stdout, setStdout] = useState('Welcome to Cis-QL Studio.');
  const [results, setResults] = useState<Record<string, any[]>>({});
  const [isRunning, setIsRunning] = useState(false);
  const [isSaving, setIsSaving] = useState(false);
  const [activeVisTab, setActiveVisTab] = useState<'track' | 'sequence'>('track');
  const [menuOpen, setMenuOpen] = useState<string | null>(null);
  const [isDarkMode, setIsDarkMode] = useState(false);
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
    const mediaQuery = window.matchMedia('(prefers-color-scheme: dark)');
    setIsDarkMode(mediaQuery.matches);

    const handler = (e: MediaQueryListEvent) => setIsDarkMode(e.matches);
    mediaQuery.addEventListener('change', handler);
    return () => mediaQuery.removeEventListener('change', handler);
  }, []);

  const defineThemes = (monacoInstance: any) => {
    monacoInstance.editor.defineTheme('cql-dark', {
      base: 'vs-dark',
      inherit: true,
      rules: [
        { token: 'keyword', foreground: '569cd6', fontStyle: 'bold' },
        { token: 'type.identifier', foreground: '4ec9b0' },
        { token: 'string', foreground: 'ce9178' },
        { token: 'number', foreground: 'b5cea8' },
        { token: 'comment', foreground: '6a9955' },
      ],
      colors: {
        'editor.background': '#1e1e1e',
        'editor.lineHighlightBackground': '#2a2d2e',
        'editorLineNumber.foreground': '#858585',
      }
    });

    monacoInstance.editor.defineTheme('cql-light', {
      base: 'vs',
      inherit: true,
      rules: [
        { token: 'keyword', foreground: '0000ff', fontStyle: 'bold' },
        { token: 'type.identifier', foreground: '267f99' },
        { token: 'string', foreground: 'a31515' },
        { token: 'number', foreground: '098658' },
        { token: 'comment', foreground: '008000' },
      ],
      colors: {
        'editor.background': '#ffffff',
        'editor.lineHighlightBackground': '#f0f0f0',
        'editorLineNumber.foreground': '#237893',
      }
    });
  };

  useEffect(() => {
    if (monaco) {
      monaco.languages.register({ id: 'cql' });
      monaco.languages.setMonarchTokensProvider('cql', cqlLanguageDef as any);
      
      defineThemes(monaco);
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
      
      folderPath = prompt("Enter the absolute path to your Cis-QL workspace directory:");
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

  return (
    <div className="app-container" onClick={() => setMenuOpen(null)}>
      {}
      <header className="top-menu-bar">
        <div className="menu-brand">Cis-QL Studio</div>
        <div className="menu-items">
          <div className="menu-item-wrapper" style={{ position: 'relative' }}>
            <div className="menu-item" onClick={(e) => { e.stopPropagation(); setMenuOpen(menuOpen === 'file' ? null : 'file'); }}>File</div>
            {menuOpen === 'file' && (
              <div className="dropdown-menu">
                <div className="dropdown-item" onClick={handleOpenFolder}>Open Folder...</div>
                <div className="dropdown-item" onClick={() => { handleSave(); setMenuOpen(null); }}>Save</div>
              </div>
            )}
          </div>
          <div className="menu-item-wrapper" style={{ position: 'relative' }}>
            <div className="menu-item" onClick={(e) => { e.stopPropagation(); setMenuOpen(menuOpen === 'run' ? null : 'run'); }}>Run</div>
            {menuOpen === 'run' && (
              <div className="dropdown-menu">
                <div className="dropdown-item" onClick={() => { handleExecute(); setMenuOpen(null); }}>Run Code</div>
              </div>
            )}
          </div>
          <div className="menu-item-wrapper" style={{ position: 'relative' }}>
            <div className="menu-item" onClick={(e) => { e.stopPropagation(); setMenuOpen(menuOpen === 'terminal' ? null : 'terminal'); }}>Terminal</div>
            {menuOpen === 'terminal' && (
              <div className="dropdown-menu">
                <div className="dropdown-item" onClick={() => { setStdout(''); setMenuOpen(null); }}>Clear Terminal</div>
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
      </header>
      
      <Split className="workspace" sizes={[20, 80]} minSize={150} gutterSize={4}>
        {}
        <div className="sidebar">
          <FileExplorer 
            onFileSelect={onFileSelect} 
            workspacePath={workspacePath} 
            onSetWorkspace={handleOpenFolder} 
          />
        </div>

        {}
        <Split className="main-content" direction="vertical" sizes={[70, 30]} minSize={100} gutterSize={4}>
          
          <Split className="editor-visualizer-split" direction="horizontal" sizes={[55, 45]} minSize={200} gutterSize={4}>
            {}
            <div className="pane editor-pane">
              <div className="pane-header">
                {activeFile ? <div className="file-tab">{activeFile.split('/').pop()}</div> : null}
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
                      fontSize: 14,
                      fontFamily: "'JetBrains Mono', 'Fira Code', monospace",
                      padding: { top: 16 },
                      scrollBeyondLastLine: false,
                      wordWrap: 'on'
                    }}
                  />
                ) : (
                  <div className="welcome-screen">
                    <h1>Cis-QL Studio</h1>
                    <div className="welcome-shortcuts">
                      <div className="shortcut"><span>Open Folder</span><span>From File Menu</span></div>
                      <div className="shortcut"><span>Save File</span><span>Cmd + S</span></div>
                      <div className="shortcut"><span>Run Code</span><span>Cmd + Enter</span></div>
                    </div>
                    {!workspacePath && (
                      <button className="run-button" style={{ marginTop: '16px' }} onClick={handleOpenFolder}>
                        Open Folder
                      </button>
                    )}
                  </div>
                )}
              </div>
            </div>

            {}
            <div className="pane visualizer-pane">
              <div className="pane-header" style={{ gap: '0' }}>
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
              <div className="pane-content">
                {activeVisTab === 'track' ? (
                  <TrackViewer results={results} />
                ) : (
                  <SequenceViewer results={results} />
                )}
              </div>
            </div>
          </Split>

          {}
          <div className="pane terminal-pane">
            <div className="pane-header">TERMINAL</div>
            <pre className="terminal-output">{stdout}</pre>
          </div>
          
        </Split>
      </Split>
    </div>
  );
}

export default App;
