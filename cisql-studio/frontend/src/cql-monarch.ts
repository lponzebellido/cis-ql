export const cqlLanguageDef = {
  keywords: [
    'LOAD', 'SEQUENCE', 'ANNOTATION', 'MATRIX', 'AS', 'FIND', 'MOTIF',
    'WITHIN', 'FROM', 'GENE', 'PROMOTER', 'ENHANCER', 'EXON', 'INTRON',
    'UTR', 'TSS', 'CDS', 'REGION', 'STRAND', 'POSITIVE', 'NEGATIVE',
    'CHR', 'EXTRACT', 'WHERE', 'INTERSECT', 'UNION', 'EXCEPT', 'AND', 'OR', 'NOT',
    'LENGTH', 'SIMILARITY', 'GC_CONTENT', 'CPG_ISLANDS', 'BP', 'KB', 'MB',
    'UPSTREAM', 'DOWNSTREAM', 'SCAN', 'THRESHOLD', 'ANALYZE', 'WINDOW',
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

      [/\d*\.\d+([eE][\-+]?\d+)?/, 'number.float'],
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
