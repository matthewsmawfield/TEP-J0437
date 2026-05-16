#!/usr/bin/env node
const chokidar = require('chokidar');
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');
const { buildStaticSite } = require('./build.js');
const { HTMLToMarkdownConverter } = require('./html-to-markdown.js');

class DevServer {
    constructor() {
        this.isBuilding = false;
        this.buildQueue = false;
        this.liveServerProcess = null;
        this.watcherReady = false;
        this.port = 51737; // Unique port for TEP-J0437
    }

    async startLiveServer() {
        console.log('🚀 Starting live server...');
        if (this.liveServerProcess) this.liveServerProcess.kill();
        this.liveServerProcess = spawn('npx', [
            'live-server', 'dist', `--port=${this.port}`, '--host=localhost', '--no-browser', '--wait=500'
        ], { stdio: 'pipe', cwd: __dirname });
        this.liveServerProcess.stdout.on('data', d => console.log('📡 ' + d.toString().trim()));
    }

    async build() {
        if (this.isBuilding) { this.buildQueue = true; return; }
        this.isBuilding = true;
        console.log('\n🔄 Rebuilding site...');
        try {
            await buildStaticSite();
            console.log('📝 Generating markdown...');
            const converter = new HTMLToMarkdownConverter();
            await converter.convertSiteToMarkdown();
            console.log('✅ Build & Markdown complete!');
            if (this.buildQueue) {
                this.buildQueue = false;
                setTimeout(() => { this.isBuilding = false; this.build(); }, 100);
                return;
            }
        } catch (e) { console.error('❌ Build failed:', e.message); }
        this.isBuilding = false;
    }

    async start() {
        console.log(`\n🎯 TEP-J0437 Development Server`);
        const distDir = path.join(__dirname, 'dist');
        if (!fs.existsSync(distDir)) fs.mkdirSync(distDir, { recursive: true });
        await this.build();
        await this.startLiveServer();
        
        const watcher = chokidar.watch([
            path.join(__dirname, 'components'),
            path.join(__dirname, 'index.html'),
            path.join(__dirname, 'manifest.json'),
            path.join(__dirname, 'figures')
        ], { ignored: ['dist/**'], persistent: true, ignoreInitial: true });

        watcher.on('change', () => this.build());
        watcher.on('add', () => this.build());
        console.log(`\n🌐 Server running at: http://localhost:${this.port}`);
    }
}

if (require.main === module) { const s = new DevServer(); s.start(); }
module.exports = DevServer;
