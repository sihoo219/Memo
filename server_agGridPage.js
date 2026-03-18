const https = require('https');
const fs = require('fs');
const express = require('express');
const sqlite3 = require('sqlite3').verbose();
const cors = require('cors');
const { default: axios } = require('axios');

const app = express();
const db = new sqlite3.Database('./gwas.db');
app.use(cors());
app.use(express.json());

const sslOptions = {
  ca: fs.readFileSync('/web_service/apache-2.4.57/conf/certs/rootca_donga.ac.kr.pem'),
  key: fs.readFileSync('/web_service/apache-2.4.57/conf/certs/prv_donga.ac.kr2.pem'),
  cert: fs.readFileSync('/web_service/apache-2.4.57/conf/certs/cert_donga.ac.kr.pem'),
}

const HTTPS_PORT = 3001; 

https.createServer(sslOptions, app).listen(HTTPS_PORT, () => {
  console.log(`node.js is running in ${HTTPS_PORT}.`);
});

// Data calling from gwas database - snps table
app.get('/api/data', (req, res) => {
  db.all("SELECT * FROM snps", (err, rows) => {
    if (err) {
      res.status(500).json({ error: err.message })
      return
    }
    const formattedRows = rows.map(row => ({
      ...row,
      p: row.p === 0 ? 1e-300 : row.p,
      logp: !isFinite(row.logp) ? 300 : Number(row.logp.toFixed(3)),
      // p: typeof row.p === 'number' ? Number(row.p.toFixed(3)) : row.p,
      // logp: typeof row.logp === 'number' ? Number(row.logp.toFixed(3)) : row.logp,
      effect: typeof row.effect === 'number' ? Number(row.effect.toFixed(3)) : row.effect
    }));
    res.json(formattedRows)
  });
});

// Data calling from samples table
app.get('/api/sample', (req, res) => {
  const trait = req.query.trait;
  console.log('Server say, trait = ', trait)
  const rows = db.all(
    'SELECT * FROM samples WHERE trait = ?', [trait], (err, rows) => {
      if (err) {
        console.error(err.message);
        res.status(500).json({ error: err.message });
        return;
      }
      // console.log('rows = ', rows);
      res.json(rows);
    }
  );
});

// Apply filtering information
app.post('/api/data/filter', (req, res) => {
  const { trait, chr, startPos, endPos, gene, pval, logp, effect } = req.body;

  let query = "SELECT * FROM snps WHERE 1=1";
  const params = [];

  if (trait) {
    query += " AND trait = ?";
    params.push(trait);
  }

  if (chr) {
    query += " AND chr = ?";
    params.push(chr);
  }

  if (startPos && endPos) {
    query += " AND pos BETWEEN ? AND ?";
    params.push(startPos, endPos)
  }

  if (gene) {
    query += " AND gene LIKE ?";
    params.push(`%${gene}%`);
  }

  if (pval !== '') {
    query += " AND p <= ?";
    params.push(parseFloat(pval));
  }

  if (logp !== '') {
    query += " AND logp >= ?";
    params.push(parseFloat(logp));
  }

  if (effect) {
    query += " AND ABS(effect) >= ?";
    params.push(parseFloat(effect));
  }

  db.all(query, params, (err, rows) => {
    if (err) {
      return res.status(500).json({ error: err.message });
    }

    const formattedRows = rows.map(row => ({
      ...row,
      // p: typeof row.p === 'number' ? Number(row.p.toFixed(3)) : row.p,
      // logp: typeof row.logp === 'number' ? Number(row.logp.toFixed(3)) : row.logp,
      p: row.p === 0 ? 1e-300 : row.p,
      logp: !isFinite(row.logp) ? 300 : Number(row.logp.toFixed(3)),
      effect: typeof row.effect === 'number' ? Number(row.effect.toFixed(3)) : row.effect,

    }));

    res.json(formattedRows);
  });
});


// KASP primer design
const { spawn } = require('child_process');

app.post('/api/marker', (req, res) => {
  const selectedData = req.body.selectedData || [];
  const py = spawn('python3', ['./Python/kasp_server_dict.py', JSON.stringify(selectedData)])
  
  let output = '';
  let errorOutput = '';

  py.stdout.on('data', (data) => {
    output += data.toString();
  });

  py.stderr.on('data', (data) => {
    errorOutput += data.toString();
  });

  py.on('close', (code) => {
    if (code !== 0) {
      console.error('Python script error:', errorOutput);
      return res.status(500).json({ error: 'python script failed', details: errorOutput });
    }
    try {
      const result = JSON.parse(output);  // Python 스크립트가 JSON 문자열 출력 가정
      res.json({ status: 'ok', result });
    } catch (err) {
      console.error('Failed to parse Python output:', err);
      res.status(500).json({ error: 'Invalid Python output' });
    }
  });
});

// app.listen(3001, () => {
//   console.log('API server running on http://10.9.0.36:3001')
// })

// Guide RNA design
app.post('/api/guide', (req, res) => {
  const selectedParams = req.body.selectedParams || [];
  const py = spawn('/home/sihoo219/anaconda3/bin/python', ['./Python/rgen_server_dict.py', JSON.stringify(selectedParams)])

  let output = '';
  let errorOutput = '';

  py.stdout.on('data', (data) => {
    output += data.toString();
  })

  py.stderr.on('data', (data) => {
    errorOutput += data.toString();
  });

  py.on('close', (code) => {
    if (code !== 0) {
      console.error('python script error: ', errorOutput);
      return res.status(500).json({ error: 'python script failed', details: errorOutput });
    }
    try {
      const result = JSON.parse(output);
      res.json({ status: 'ok', result });
    } catch (err) {
      console.error('Failed to parse Python output:', err);
      res.status(500).json({ error: 'Invalid Python output' });
    }
  })
});

// app.listen(3001, () => {
//   console.log('API server running on http://10.9.0.36:3001')
// })

// Find near snps
app.post('/api/near-snps', (req, res) => {
  const selectedInfo = req.body.info
  console.log('Server say: ', selectedInfo)
  const py = spawn('python3', ['./Python/nearSnpGeneAnn_server_dict.py', JSON.stringify(selectedInfo)])

  let output = '';
  let errorOutput = '';

  py.stdout.on('data', (data) => {
    output += data.toString();
  });

  py.stderr.on('data', (data) => {
    errorOutput += data.toString();
  });

  py.on('close', (code) => {
    if (code !== 0) {
      console.error('Python script error:', errorOutput);
      return res.status(500).json({ error: 'python script failed', details: errorOutput });
    }
    try {
      const result = JSON.parse(output);
      res.json({ status: 'ok', result });
    } catch (err) {
      console.error('Failed to parse Python output:', err);
      res.status(500).json({ error: 'Invalid Python output' });
    }
  });

})

// app.listen(3001, () => {
//   console.log('API server running on http://10.9.0.36:3001')
// })


// Haplotype (continious)
app.post('/api/haplotypes', (req, res) => {
  const selectedInfo = req.body
  // console.log('Server say: ', selectedInfo.type)
  const py = spawn('python3', ['./Python/haplotype_server_dict.py', JSON.stringify(selectedInfo)])

  let output = '';
  let errorOutput = '';

  py.stdout.on('data', (data) => {
    output += data.toString();
  });

  py.stderr.on('data', (data) => {
    errorOutput += data.toString();
  });

  py.on('close', (code) => {
    if (code !== 0) {
      console.error('Python script error:', errorOutput);
      return res.status(500).json({ error: 'python script failed', details: errorOutput });
    }
    try {
      const result = JSON.parse(output);
      res.json({ status: 'ok', result });
    } catch (err) {
      console.error('Failed to parse Python output:', err);
      res.status(500).json({ error: 'Invalid Python output' });
    }
  });

})

// app.listen(3001, () => {
//   console.log('API server running on http://10.9.0.36:3001')
// })