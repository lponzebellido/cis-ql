import React, { useState, useEffect, useCallback } from 'react';

interface FileNode {
  name: string;
  isDirectory: boolean;
  path: string;
  children?: FileNode[];
}

interface FileExplorerProps {
  onFileSelect: (path: string, content: string) => void;
  workspacePath: string | null;
  onSetWorkspace: () => void;
}

const EXT_COLORS: Record<string, string> = {
  cql: '#3fb950',
  fasta: '#58a6ff',
  fa: '#58a6ff',
  gff3: '#d29922',
  gff: '#d29922',
  pwm: '#bc8cff',
  txt: '#8b949e',
  md: '#8b949e',
  cpp: '#f85149',
  h: '#f85149',
  js: '#e3b341',
  ts: '#58a6ff',
  tsx: '#58a6ff',
  json: '#d29922',
};

function getFileColor(name: string): string {
  const ext = name.split('.').pop()?.toLowerCase() || '';
  return EXT_COLORS[ext] || '#484f58';
}

export const FileExplorer: React.FC<FileExplorerProps> = ({ onFileSelect, workspacePath, onSetWorkspace }) => {
  const [nodes, setNodes] = useState<FileNode[]>([]);
  const [openPaths, setOpenPaths] = useState<Set<string>>(new Set());
  const [selectedPath, setSelectedPath] = useState<string | null>(null);
  
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  
  const [creatingType, setCreatingType] = useState<'file' | 'dir' | null>(null);
  const [creatingParent, setCreatingParent] = useState('');
  const [creatingName, setCreatingName] = useState('');
  
  const [renamingPath, setRenamingPath] = useState<string | null>(null);
  const [renamingName, setRenamingName] = useState('');

  const [contextMenu, setContextMenu] = useState<{x: number, y: number, node: FileNode} | null>(null);

  const fetchDir = async (path: string = '') => {
    try {
      const res = await fetch(`http://localhost:3001/api/fs/list?path=${encodeURIComponent(path)}`);
      if (!res.ok) throw new Error('Failed to fetch directory');
      return await res.json() as FileNode[];
    } catch (err: any) {
      console.error(err);
      return [];
    }
  };

  const loadTree = useCallback(async () => {
    if (!workspacePath) return;
    
    async function populateChildren(nodeList: FileNode[]): Promise<FileNode[]> {
      const result: FileNode[] = [];
      for (const node of nodeList) {
        let children = node.children;
        if (node.isDirectory && openPaths.has(node.path)) {
          const freshChildren = await fetchDir(node.path);
          children = await populateChildren(freshChildren);
        }
        result.push({ ...node, children });
      }
      return result;
    }

    try {
      const rootNodes = await fetchDir('');
      const populatedNodes = await populateChildren(rootNodes);
      setNodes(populatedNodes);
      setError(null);
    } catch (err: any) {
      setError(err.message);
    }
  }, [workspacePath, openPaths]);

  useEffect(() => {
    if (workspacePath) {
      setLoading(true);
      loadTree().then(() => setLoading(false));
    } else {
      setNodes([]);
      setOpenPaths(prev => prev.size === 0 ? prev : new Set());
      setSelectedPath(null);
    }
  }, [workspacePath, loadTree]);

  useEffect(() => {
    const closeMenu = () => setContextMenu(null);
    window.addEventListener('click', closeMenu);
    return () => window.removeEventListener('click', closeMenu);
  }, []);

  const handleToggle = async (node: FileNode) => {
    setSelectedPath(node.path);

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

    setOpenPaths(prev => {
      const newPaths = new Set(prev);
      if (newPaths.has(node.path)) {
        newPaths.delete(node.path);
      } else {
        newPaths.add(node.path);
      }
      return newPaths;
    });
  };

  const findNode = (nodeList: FileNode[], path: string): FileNode | null => {
    for (const node of nodeList) {
      if (node.path === path) return node;
      if (node.children) {
        const found = findNode(node.children, path);
        if (found) return found;
      }
    }
    return null;
  };

  const getParentPath = (pathStr: string) => {
    const parts = pathStr.split('/');
    parts.pop();
    return parts.join('/');
  };

  const startCreating = (type: 'file' | 'dir') => {
    let targetParent = '';
    
    if (selectedPath) {
      const selNode = findNode(nodes, selectedPath);
      if (selNode && selNode.isDirectory) {
        targetParent = selNode.path;
        setOpenPaths(prev => {
          const np = new Set(prev);
          np.add(targetParent);
          return np;
        });
      } else if (selNode && !selNode.isDirectory) {
        targetParent = getParentPath(selNode.path);
      }
    }

    setCreatingType(type);
    setCreatingParent(targetParent);
    setCreatingName('');
  };

  const submitCreate = async () => {
    if (!creatingType) return;
    const type = creatingType;
    const name = creatingName;
    const parentPath = creatingParent;
    setCreatingType(null);
    if (!name.trim()) return;

    const targetPath = parentPath ? `${parentPath}/${name}` : name;
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
      
      setOpenPaths(prev => {
        const np = new Set(prev);
        if (parentPath) np.add(parentPath);
        if (type === 'dir') np.add(targetPath);
        return np;
      });
      setSelectedPath(targetPath);
      await loadTree();
      
      if (type === 'file') {
        onFileSelect(targetPath, '');
      }
    } catch (err: any) {
      alert(err.message);
    }
  };

  const handleDragStart = (e: React.DragEvent, path: string) => {
    e.dataTransfer.setData('application/x-cisql-filepath', path);
    e.dataTransfer.effectAllowed = 'move';
  };

  const handleDragOver = (e: React.DragEvent, node: FileNode) => {
    if (node.isDirectory) {
      e.preventDefault();
      e.dataTransfer.dropEffect = 'move';
    }
  };

  const handleDrop = async (e: React.DragEvent, targetNode: FileNode) => {
    e.preventDefault();
    if (!targetNode.isDirectory) return;

    const oldPath = e.dataTransfer.getData('application/x-cisql-filepath');
    if (!oldPath || oldPath === targetNode.path) return;

    const fileName = oldPath.split('/').pop();
    const targetParentPath = targetNode.path === '' ? '' : targetNode.path + '/';
    const newPath = targetParentPath + fileName;

    if (oldPath === newPath) return;

    try {
      const res = await fetch('http://localhost:3001/api/fs/rename', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ oldPath, newPath })
      });
      if (!res.ok) throw new Error((await res.json()).error);
      
      setOpenPaths(prev => {
        const np = new Set(prev);
        np.add(targetNode.path);
        return np;
      });
      
      if (selectedPath === oldPath) setSelectedPath(newPath);
      await loadTree();
    } catch (err: any) {
      alert('Failed to move: ' + err.message);
    }
  };

  const handleContextMenu = (e: React.MouseEvent, node: FileNode) => {
    e.preventDefault();
    setSelectedPath(node.path);
    setContextMenu({ x: e.clientX, y: e.clientY, node });
  };

  const startRenaming = (node: FileNode) => {
    setRenamingPath(node.path);
    setRenamingName(node.name);
  };

  const submitRename = async () => {
    if (!renamingPath || !renamingName.trim()) {
      setRenamingPath(null);
      return;
    }
    const node = findNode(nodes, renamingPath);
    if (!node || node.name === renamingName) {
      setRenamingPath(null);
      return;
    }

    const parentPath = getParentPath(renamingPath);
    const newPath = parentPath ? `${parentPath}/${renamingName}` : renamingName;

    try {
      const res = await fetch('http://localhost:3001/api/fs/rename', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ oldPath: renamingPath, newPath })
      });
      if (!res.ok) throw new Error((await res.json()).error);
      
      setOpenPaths(prev => {
        const np = new Set(prev);
        if (np.has(renamingPath)) {
          np.delete(renamingPath);
          np.add(newPath);
        }
        return np;
      });
      
      if (selectedPath === renamingPath) setSelectedPath(newPath);
      setRenamingPath(null);
      await loadTree();
    } catch (err: any) {
      alert('Rename failed: ' + err.message);
      setRenamingPath(null);
    }
  };

  const handleDelete = async (node: FileNode) => {
    if (!confirm(`Are you sure you want to delete '${node.name}'?`)) return;
    try {
      const res = await fetch('http://localhost:3001/api/fs/delete', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ path: node.path })
      });
      if (!res.ok) throw new Error((await res.json()).error);
      
      if (selectedPath === node.path) setSelectedPath(null);
      await loadTree();
    } catch (err: any) {
      alert('Delete failed: ' + err.message);
    }
  };

  const renderNodes = (nodesToRender: FileNode[], depth: number = 0) => {
    return nodesToRender.map(node => {
      const isOpen = openPaths.has(node.path);
      const isSelected = selectedPath === node.path;
      const isRenaming = renamingPath === node.path;

      return (
        <div key={node.path} className="file-node-container">
          <div
            className={`file-node ${node.isDirectory ? 'directory' : 'file'} ${isSelected ? 'active' : ''}`}
            onClick={(e) => {
              e.stopPropagation();
              if (!isRenaming) handleToggle(node);
            }}
            onContextMenu={(e) => handleContextMenu(e, node)}
            draggable={!isRenaming}
            onDragStart={(e) => handleDragStart(e, node.path)}
            onDragOver={(e) => handleDragOver(e, node)}
            onDrop={(e) => handleDrop(e, node)}
            style={{ paddingLeft: `${(depth + 1) * 12}px` }}
          >
            {node.isDirectory ? (
              <span className="file-icon" style={{ color: 'var(--text-muted)' }}>
                {isOpen ? '▾' : '▸'}
              </span>
            ) : (
              <span className="file-dot" style={{ background: getFileColor(node.name) }}></span>
            )}
            
            {isRenaming ? (
              <input
                className="rename-input"
                autoFocus
                value={renamingName}
                onChange={e => setRenamingName(e.target.value)}
                onBlur={submitRename}
                onKeyDown={e => {
                  if (e.key === 'Enter') submitRename();
                  if (e.key === 'Escape') setRenamingPath(null);
                }}
                onClick={e => e.stopPropagation()}
              />
            ) : (
              <span>{node.name}</span>
            )}
          </div>
          
          {node.isDirectory && isOpen && (
            <div className="file-children" style={{ marginLeft: `${(depth + 1) * 12 + 4}px` }}>
              {creatingType && creatingParent === node.path && (
                <div className="new-item-input">
                  <input
                    type="text"
                    autoFocus
                    placeholder={creatingType === 'dir' ? 'folder name' : 'file name'}
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
              {node.children && renderNodes(node.children, depth + 1)}
            </div>
          )}
        </div>
      );
    });
  };

  if (!workspacePath) {
    return (
      <div className="file-explorer" style={{ padding: '24px 16px', textAlign: 'center' }}>
        <p style={{ color: 'var(--text-muted)', fontSize: '0.8rem', marginBottom: '16px' }}>No folder opened</p>
        <button className="welcome-open-btn" onClick={onSetWorkspace}>
          Open Folder
        </button>
      </div>
    );
  }

  return (
    <div 
      className="file-explorer" 
      onClick={() => setSelectedPath(null)}
      onDragOver={(e) => {
        e.preventDefault();
        e.dataTransfer.dropEffect = 'move';
      }}
      onDrop={(e) => {
        e.preventDefault();
        const oldPath = e.dataTransfer.getData('application/x-cisql-filepath');
        if (!oldPath || oldPath.indexOf('/') === -1) return; 
        
        const fileName = oldPath.split('/').pop();
        const newPath = fileName as string;
        
        if (oldPath === newPath) return;

        fetch('http://localhost:3001/api/fs/rename', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ oldPath, newPath })
        }).then(res => {
          if (!res.ok) throw new Error('Failed to move');
          loadTree();
        }).catch(err => alert(err.message));
      }}
    >
      <div className="explorer-header">
        <span>{workspacePath.split('/').pop()?.toUpperCase()}</span>
        <div className="explorer-actions" onClick={e => e.stopPropagation()}>
          <button onClick={() => startCreating('file')} title="New File">+</button>
          <button onClick={() => startCreating('dir')} title="New Folder">&#x2750;</button>
        </div>
      </div>
      <div className="explorer-content">
        {loading && nodes.length === 0 ? (
          <div style={{ padding: '12px 16px', color: 'var(--text-muted)', fontSize: '0.8rem' }}>Loading...</div>
        ) : error ? (
          <div style={{ padding: '12px 16px', color: 'var(--red)', fontSize: '0.8rem' }}>{error}</div>
        ) : (
          <>
            {creatingType && creatingParent === '' && (
              <div className="new-item-input">
                <input
                  type="text"
                  autoFocus
                  placeholder={creatingType === 'dir' ? 'folder name' : 'file name'}
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

      {contextMenu && (
        <div 
          className="context-menu dropdown-menu" 
          style={{ top: contextMenu.y, left: contextMenu.x, position: 'fixed' }}
          onClick={e => e.stopPropagation()}
        >
          <div className="dropdown-item" onClick={() => { startRenaming(contextMenu.node); setContextMenu(null); }}>
            Rename
          </div>
          <div className="dropdown-item" style={{ color: 'var(--red)' }} onClick={() => { handleDelete(contextMenu.node); setContextMenu(null); }}>
            Delete
          </div>
        </div>
      )}
    </div>
  );
};
