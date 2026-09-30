type AnnotationEvidence = {
  source: string;
  score: string;
  phase: string;
  id: string;
  name: string;
  parents: string[];
  attributes: Record<string, string>;
};

type FeatureGroupMember = {
  index: number;
  chr: string;
  start: number;
  end: number;
  length: number;
  strand: string;
  type: string;
  name: string;
  annotationEvidence?: AnnotationEvidence;
};

type FeatureGroup = {
  sourceAlias: string;
  id: string;
  chr: string;
  type: string;
  strand: string;
  span: {
    start: number;
    end: number;
    length: number;
  };
  totalMemberLength: number;
  memberCount: number;
  members: FeatureGroupMember[];
};

type FeatureGroupViewerProps = {
  collections: Record<string, FeatureGroup[]>;
};

export function FeatureGroupViewer({ collections }: FeatureGroupViewerProps) {
  const entries = Object.entries(collections);
  if (entries.length === 0) {
    return (
      <div className="feature-group-empty">
        Run <code>GROUP annotation BY ID AS groups;</code> to inspect logical GFF3 features.
      </div>
    );
  }

  return (
    <div className="feature-group-viewer">
      {entries.map(([name, groups]) => (
        <section className="feature-group-collection" key={name}>
          <div className="feature-group-collection-heading">
            <div>
              <div className="feature-group-collection-name">{name}</div>
              <div className="feature-group-collection-count">{groups.length} groups</div>
            </div>
          </div>
          {groups.map((group) => (
            <article className="feature-group" key={group.id}>
              <div className="feature-group-heading">
                <div>
                  <div className="feature-group-id">{group.id}</div>
                  <div className="feature-group-location">
                    {group.chr}:{group.span.start}-{group.span.end} · {group.type} · {group.strand}
                  </div>
                </div>
                <span>{group.sourceAlias}</span>
              </div>
              <div className="feature-group-metrics">
                <div><span>Members</span><strong>{group.memberCount}</strong></div>
                <div><span>Member length</span><strong>{group.totalMemberLength}</strong></div>
                <div><span>Genomic span</span><strong>{group.span.length}</strong></div>
                <div><span>Span − members</span><strong>{group.span.length - group.totalMemberLength}</strong></div>
              </div>
              <table>
                <thead>
                  <tr><th>#</th><th>Start</th><th>End</th><th>Length</th><th>Score</th><th>Phase</th><th>Parents</th></tr>
                </thead>
                <tbody>
                  {group.members.map((member) => (
                    <tr key={`${group.id}-${member.index}`}>
                      <td>{member.index}</td>
                      <td>{member.start}</td>
                      <td>{member.end}</td>
                      <td>{member.length}</td>
                      <td>{member.annotationEvidence?.score || ''}</td>
                      <td>{member.annotationEvidence?.phase || ''}</td>
                      <td>{member.annotationEvidence?.parents.join(', ') || ''}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </article>
          ))}
        </section>
      ))}
    </div>
  );
}
