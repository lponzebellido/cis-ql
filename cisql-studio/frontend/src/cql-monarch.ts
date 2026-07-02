export const cqlLanguageDef = {
  // Set defaultToken to invalid to see what you do not tokenize yet
  // defaultToken: 'invalid',

  keywords: [
    'LOAD', 'SEQUENCE', 'ANNOTATION', 'MATRIX', 'AS', 'FIND', 'MOTIF',
    'WITHIN', 'FROM', 'GENE', 'PROMOTER', 'ENHANCER', 'EXON', 'INTRON',
    'UTR', 'TSS', 'CDS', 'REGION', 'STRAND', 'POSITIVE', 'NEGATIVE',
    'CHR', 'EXTRACT', 'WHERE', 'INTERSECT', 'UNION', 'EXCEPT', 'AND', 'OR', 'NOT',
    'LENGTH', 'SIMILARITY', 'BP', 'KB', 'MB', 'UPSTREAM', 'DOWNSTREAM', 'SCAN', 'THRESHOLD'
  ],

  operators: [
    '<=', '>=', '==', '!=', '>', '<', '=', '%'
  ],

  // we include these common regular expressions
  symbols:  /[=><!~?:&|+\-*\/\^%]+/,

  // The main tokenizer for our languages
  tokenizer: {
    root: [
      // identifiers and keywords
      [/[a-z_$][\w$]*/, { cases: { '@keywords': 'keyword',
                                   '@default': 'identifier' } }],
      [/[A-Z][\w\$]*/, { cases: { '@keywords': 'keyword',
                                  '@default': 'type.identifier' } }],

      // whitespace
      { include: '@whitespace' },

      // delimiters and operators
      [/[{}()\[\]]/, '@brackets'],
      [/[<>](?!@symbols)/, '@brackets'],
      [/@symbols/, { cases: { '@operators': 'operator',
                              '@default'  : '' } } ],

      // numbers
      [/\d*\.\d+([eE][\-+]?\d+)?/, 'number.float'],
      [/0[xX][0-9a-fA-F]+/, 'number.hex'],
      [/\d+/, 'number'],

      // delimiter: after number because of .\d floats
      [/[;,.]/, 'delimiter'],

      // strings
      [/"([^"\\]|\\.)*$/, 'string.invalid' ],  // non-teminated string
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
      [/\/\*/,    'comment', '@push' ],    // nested comment
      ["\\*/",    'comment', '@pop'  ],
      [/[\/*]/,   'comment' ]
    ],
  },
};
