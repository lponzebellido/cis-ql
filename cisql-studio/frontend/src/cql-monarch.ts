export const cqlLanguageDef = {
  keywords: [
    'LOAD', 'USE', 'EXPORT', 'DEFINE', 'SEQUENCE', 'ANNOTATION', 'TRACK', 'MATRIX',
    'AS', 'OF', 'PROMOTERS', 'FIND', 'MOTIF', 'MODULE',
    'WITHIN', 'WITH', 'FROM', 'TO', 'GENE', 'PROMOTER', 'ENHANCER', 'EXON', 'INTRON',
    'UTR', 'TSS', 'CDS', 'REGION', 'STRAND', 'POSITIVE', 'NEGATIVE',
    'CHR', 'EXTRACT', 'WHERE', 'INTERSECT', 'UNION', 'EXCEPT', 'OVERLAPS', 'NEAR',
    'CONSENSUS', 'ANCHOR', 'MIN_SUPPORT', 'SUPPORT_COUNT',
    'AND', 'OR', 'NOT', 'COUNT', 'SPACING', 'ORDER', 'ORIENTATION',
    'ANY', 'SAME', 'OPPOSITE', 'AS_WRITTEN',
    'EVIDENCE', 'ACCESSIBILITY', 'BINDING', 'OTHER', 'ASSAY', 'SAMPLE',
    'CONDITION', 'REPLICATE', 'CONTROL',
    'TRACK_SCORE', 'SIGNAL_VALUE', 'MINUS_LOG10_PVALUE',
    'MINUS_LOG10_QVALUE', 'EVIDENCE_CLASS',
    'LENGTH', 'SIMILARITY', 'GC_CONTENT', 'CPG_ISLANDS', 'BP', 'KB', 'MB',
    'UPSTREAM', 'DOWNSTREAM', 'SCAN', 'THRESHOLD', 'PVALUE', 'QVALUE',
    'BACKGROUND', 'UNIFORM', 'ANALYZE', 'WINDOW', 'FORMAT', 'BED', 'NARROWPEAK', 'GFF3', 'TSV',
    'IF', 'THEN', 'ELSE', 'ENDIF', 'FOREACH', 'IN', 'DO', 'ENDFOR'
  ],

  operators: [
    '<=', '>=', '>', '<', '=', '%'
  ],

  symbols:  /[=><!~?:&|+\-*\/\^%]+/,

  tokenizer: {
    root: [
      [/[a-z_$][\w$]*/, { cases: { '@keywords': 'keyword',
                                   '@default': 'identifier' } }],
      [/[A-Z][\w\$]*/, { cases: { '@keywords': 'keyword',
                                  '@default': 'type.identifier' } }],

      { include: '@whitespace' },

      [/[{}()\[\]]/, '@brackets'],
      [/[<>](?!@symbols)/, '@brackets'],
      [/@symbols/, { cases: { '@operators': 'operator',
                              '@default'  : '' } } ],

      [/\d+(\.\d*)?[eE][\-+]?\d+/, 'number.float'],
      [/\d*\.\d+/, 'number.float'],
      [/0[xX][0-9a-fA-F]+/, 'number.hex'],
      [/\d+/, 'number'],

      [/[;,.]/, 'delimiter'],

      [/"([^"\\]|\\.)*$/, 'string.invalid' ],
      [/"/,  { token: 'string.quote', bracket: '@open', next: '@string' } ],
    ],

    string: [
      [/[^\\"]+/,  'string'],
      [/\\./,      'string.escape.invalid'],
      [/"/,        { token: 'string.quote', bracket: '@close', next: '@pop' } ]
    ],

    whitespace: [
      [/[ \t\r\n]+/, 'white'],
      [/\/\*/,       'comment', '@comment' ],
      [/\/\/.*$/,    'comment'],
    ],

    comment: [
      [/[^\/*]+/, 'comment' ],
      [/\/\*/,    'comment', '@push' ],
      ["\\*/",    'comment', '@pop'  ],
      [/[\/*]/,   'comment' ]
    ],
  },
};
