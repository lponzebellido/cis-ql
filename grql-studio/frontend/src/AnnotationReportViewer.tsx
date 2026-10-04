type MultiRecordIdentity = {
  id: string;
  recordNumbers: number[];
};

type UnresolvedParent = {
  recordNumber: number;
  chr: string;
  type: string;
  childId: string;
  childName: string;
  parentId: string;
};

type IdentityConflict = {
  id: string;
  recordNumbers: number[];
  chromosomes: string[];
  types: string[];
  strands: string[];
};

type AnnotationReport = {
  annotationAlias: string;
  valid: boolean;
  summary: Record<string, number>;
  multiRecordIdentities: MultiRecordIdentity[];
  unresolvedParents: UnresolvedParent[];
  identityConflicts: IdentityConflict[];
  cycles: string[][];
};

type AnnotationReportViewerProps = {
  reports: Record<string, AnnotationReport>;
};

const SUMMARY_LABELS: Record<string, string> = {
  totalRecords: 'Records',
  recordsWithId: 'Records with ID',
  uniqueIds: 'Unique IDs',
  parentReferences: 'Parent references',
  multiRecordIds: 'Multi-record IDs',
  unresolvedParents: 'Unresolved parents',
  identityConflicts: 'Identity conflicts',
  cycles: 'Cycles',
};

export function AnnotationReportViewer({ reports }: AnnotationReportViewerProps) {
  const entries = Object.entries(reports);
  if (entries.length === 0) {
    return (
      <div className="annotation-report-empty">
        Run <code>VALIDATE ANNOTATION dataset AS report;</code> to inspect a GFF3 hierarchy.
      </div>
    );
  }

  return (
    <div className="annotation-report-viewer">
      {entries.map(([name, report]) => (
        <section className="annotation-report" key={name}>
          <div className="annotation-report-heading">
            <div>
              <div className="annotation-report-name">{name}</div>
              <div className="annotation-report-source">Annotation: {report.annotationAlias}</div>
            </div>
            <span className={`annotation-report-status ${report.valid ? 'valid' : 'invalid'}`}>
              {report.valid ? 'Valid' : 'Issues found'}
            </span>
          </div>

          <div className="annotation-report-summary">
            {Object.entries(report.summary).map(([key, value]) => (
              <div className="annotation-summary-card" key={key}>
                <span>{SUMMARY_LABELS[key] || key}</span>
                <strong>{value}</strong>
              </div>
            ))}
          </div>

          {report.multiRecordIdentities.length > 0 && (
            <div className="annotation-diagnostic-group">
              <h3>Multi-record identities <span>Information</span></h3>
              <table>
                <thead><tr><th>ID</th><th>Records</th></tr></thead>
                <tbody>
                  {report.multiRecordIdentities.map((item) => (
                    <tr key={item.id}><td>{item.id}</td><td>{item.recordNumbers.join(', ')}</td></tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}

          {report.unresolvedParents.length > 0 && (
            <div className="annotation-diagnostic-group error">
              <h3>Unresolved parents <span>Error</span></h3>
              <table>
                <thead><tr><th>Record</th><th>Feature</th><th>Parent</th><th>Location</th></tr></thead>
                <tbody>
                  {report.unresolvedParents.map((item, index) => (
                    <tr key={`${item.recordNumber}-${item.parentId}-${index}`}>
                      <td>{item.recordNumber}</td>
                      <td>{item.childId || item.childName || item.type}</td>
                      <td>{item.parentId}</td>
                      <td>{item.chr}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}

          {report.identityConflicts.length > 0 && (
            <div className="annotation-diagnostic-group error">
              <h3>Identity conflicts <span>Error</span></h3>
              <table>
                <thead><tr><th>ID</th><th>Records</th><th>Chromosomes</th><th>Types</th><th>Strands</th></tr></thead>
                <tbody>
                  {report.identityConflicts.map((item) => (
                    <tr key={item.id}>
                      <td>{item.id}</td>
                      <td>{item.recordNumbers.join(', ')}</td>
                      <td>{item.chromosomes.join(', ')}</td>
                      <td>{item.types.join(', ')}</td>
                      <td>{item.strands.join(', ')}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}

          {report.cycles.length > 0 && (
            <div className="annotation-diagnostic-group error">
              <h3>Cycles <span>Error</span></h3>
              <ul className="annotation-cycle-list">
                {report.cycles.map((cycle, index) => (
                  <li key={`${cycle.join('-')}-${index}`}>{cycle.join(' → ')}</li>
                ))}
              </ul>
            </div>
          )}
        </section>
      ))}
    </div>
  );
}
