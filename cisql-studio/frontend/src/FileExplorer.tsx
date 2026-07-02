import React, { useState, useEffect } from 'react';

interface FileNode {
  name: string;
  isDirectory: boolean;
  path: string;
  children?: FileNode[];
  isOpen?: boolean;
}

interface FileExplorerProps {
  onFileSelect: (path: string, content: string) => void;
  workspacePath: string | null;
  onSetWorkspace: () => void;
}

export const FileExplorer: React.FC<FileExplorerProps> = ({ onFileSelect, workspacePath, onSetWorkspace }) => {
  const [nodes, setNodes] = useState<FileNode[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  
  const [creatingType, setCreatingType] = useState<'file' | 'dir' | null>(null);
  const [creatingName, setCreatingName] = useState('');
  const [creatingPath, setCreatingPath] = useState(''); 

  const fetchDir = async (path: string = '') => {
    try {
      const res = await fetch(`http://localhost:3001/api/fs/list?path=${encodeURIComponent(path)}`);
      if (!res.ok) throw new Error('Failed to fetch directory');
      return await res.json();
    } catch (err: any) {
      setError(err.message);
      return [];
    }
  };

  const loadWorkspace = () => {
    if (!workspacePath) return;
    setLoading(true);
    fetchDir('').then(data => {
      setNodes(data);
      setLoading(false);
    });
  };

  useEffect(() => {
    if (workspacePath) {
      loadWorkspace();
    } else {
      setNodes([]);
    }
  }, [workspacePath]);

  const handleToggle = async (node: FileNode, indexPath: number[]) => {
    if (!node.isDirectory) {
      try {
        const res = await fetch(`http://localhost:3001/api/fs/read?path=${encodeURIComponent(node.path)}`);
        if (!res.ok) throw new Error('Failed to read file');
        const { content } = await res.json();
        onFileSelect(node.path, content);
      } catch (err: any) {
        alert(err.message);
      }
      return;
    }

    const newNodes = [...nodes];
    let curr = newNodes;
    for (let i = 0; i < indexPath.length - 1; i++) {
      curr = curr[indexPath[i]].children!;
    }
    const target = curr[indexPath[indexPath.length - 1]];
    
    if (target.isOpen) {
      target.isOpen = false;
      setNodes(newNodes);
    } else {
      target.isOpen = true;
      if (!target.children) {
        const children = await fetchDir(target.path);
        target.children = children;
      }
      setNodes(newNodes);
    }
  };

  const startCreating = (type: 'file' | 'dir', parentPath: string = '') => {
    setCreatingType(type);
    setCreatingPath(parentPath);
    setCreatingName('');
  };

  const submitCreate = async () => {
    if (!creatingType) return; 
    const type = creatingType;
    const name = creatingName;
    const path = creatingPath;
    
    setCreatingType(null); 
    
    if (!name.trim()) return;

    const targetPath = path ? `${path}/${name}` : name;
    const endpoint = type === 'file' ? '/api/fs/create-file' : '/api/fs/create-dir';

    try {
      const res = await fetch(`http://localhost:3001${endpoint}`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ path: targetPath })
      });
      if (!res.ok) {
        const err = await res.json();
        throw new Error(err.error || 'Failed to create');
      }
      
      
      loadWorkspace();
      
      if (type === 'file') {
        onFileSelect(targetPath, ''); 
      }
    } catch (err: any) {
      alert(err.message);
    }
  };

  const renderNodes = (nodesToRender: FileNode[], parentPath: number[] = []) => {
    return nodesToRender.map((node, index) => {
      const currentPath = [...parentPath, index];

      return (
        <div key={node.path} className="file-node-container">
          <div 
            className={`file-node ${node.isDirectory ? 'directory' : 'file'}`}
            onClick={() => handleToggle(node, currentPath)}
          >
            <span className="file-icon" style={{ display: 'inline-block', width: '16px', textAlign: 'center' }}>
              {node.isDirectory ? (node.isOpen ? '▼' : '▶') : '≡'}
            </span>
            <span className="file-name">{node.name}</span>
          </div>
          {node.isDirectory && node.isOpen && node.children && (
            <div className="file-children">
              {creatingType && creatingPath === node.path && (
                <div className="file-node new-item-input">
                  <span className="file-icon">{creatingType === 'dir' ? '▶' : '≡'}</span>
                  <input
                    type="text"
                    autoFocus
                    value={creatingName}
                    onChange={e => setCreatingName(e.target.value)}
                    onBlur={submitCreate}
                    onKeyDown={e => {
                      if (e.key === 'Enter') submitCreate();
                      if (e.key === 'Escape') setCreatingType(null);
                    }}
                  />
                </div>
              )}
              {renderNodes(node.children, currentPath)}
            </div>
          )}
        </div>
      );
    });
  };

  if (!workspacePath) {
    return (
      <div className="file-explorer" style={{ padding: '20px', textAlign: 'center' }}>
        <p style={{ color: 'var(--text-muted)', fontSize: '0.85rem', marginBottom: '16px' }}>No folder opened</p>
        <button className="run-button" style={{ margin: '0 auto' }} onClick={onSetWorkspace}>
          Open Folder
        </button>
      </div>
    );
  }

  return (
    <div className="file-explorer">
      <div className="explorer-header">
        <span>WORKSPACE</span>
        <div className="explorer-actions">
          <button onClick={() => startCreating('file')} title="New File">+ File</button>
          <button onClick={() => startCreating('dir')} title="New Folder">+ Dir</button>
        </div>
      </div>
      <div className="explorer-content">
        {loading && nodes.length === 0 ? (
          <div className="loading">Loading...</div>
        ) : error ? (
          <div className="error">{error}</div>
        ) : (
          <>
            {creatingType && creatingPath === '' && (
              <div className="file-node new-item-input">
                <span className="file-icon">{creatingType === 'dir' ? '▶' : '≡'}</span>
                <input
                  type="text"
                  autoFocus
                  value={creatingName}
                  onChange={e => setCreatingName(e.target.value)}
                  onBlur={submitCreate}
                  onKeyDown={e => {
                    if (e.key === 'Enter') submitCreate();
                    if (e.key === 'Escape') setCreatingType(null);
                  }}
                />
              </div>
            )}
            {renderNodes(nodes)}
          </>
        )}
      </div>
    </div>
  );
};
