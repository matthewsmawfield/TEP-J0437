#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const yaml = require('js-yaml');

class HTMLToMarkdownConverter {
    constructor() { this.output = ''; }

    decodeEntities(text) {
        return text
            .replace(/&amp;/g, '&')
            .replace(/&lt;/g, '<')
            .replace(/&gt;/g, '>')
            .replace(/&nbsp;/g, ' ')
            .replace(/&mdash;/g, '—')
            .replace(/&ndash;/g, '–')
            .replace(/&#10004;/g, '✔')
            .replace(/&#10008;/g, '✘');
    }

    cleanInlineHtml(html) {
        return this.decodeEntities(
            html
                .replace(/<(strong|b)[^>]*>([\s\S]*?)<\/(strong|b)>/gi, '**$2**')
                .replace(/<(em|i)[^>]*>([\s\S]*?)<\/(em|i)>/gi, '*$2*')
                .replace(/<code[^>]*>([\s\S]*?)<\/code>/gi, '`$1`')
                .replace(/<br\s*\/?>/gi, ' ')
                .replace(/<\/?[a-zA-Z][^>]*>/g, '')
        )
            .replace(/\s+/g, ' ')
            .trim();
    }

    tableToMarkdown(tableHtml) {
        const captionMatch = tableHtml.match(/<caption[^>]*>([\s\S]*?)<\/caption>/i);
        const caption = captionMatch ? this.cleanInlineHtml(captionMatch[1]) : '';
        const rows = [];
        const rowRegex = /<tr[^>]*>([\s\S]*?)<\/tr>/gi;
        let rowMatch;

        while ((rowMatch = rowRegex.exec(tableHtml)) !== null) {
            const cells = [];
            const cellRegex = /<t[hd][^>]*>([\s\S]*?)<\/t[hd]>/gi;
            let cellMatch;

            while ((cellMatch = cellRegex.exec(rowMatch[1])) !== null) {
                cells.push(this.cleanInlineHtml(cellMatch[1]).replace(/\|/g, '\\|'));
            }

            if (cells.length) rows.push(cells);
        }

        if (!rows.length) return '';

        const header = rows[0];
        const normalizeRow = (row) => `| ${header.map((_, index) => row[index] || '').join(' | ')} |`;

        let markdown = '';
        if (caption) markdown += `\n\n${caption}\n\n`;
        markdown += `${normalizeRow(header)}\n`;
        markdown += `| ${header.map(() => '---').join(' | ')} |\n`;
        for (const row of rows.slice(1)) markdown += `${normalizeRow(row)}\n`;
        return `\n\n${markdown}\n`;
    }

    cleanParagraphWhitespace(html) {
        // Remove leading indentation from content inside <p> tags
        return html.replace(/<p(\s[^>]*)?>([\s\S]*?)<\/p>/gi, (match, attrs, content) => {
            if (content.match(/^\s/)) {
                const lines = content.split('\n');
                const cleanedLines = lines.map(line => line.replace(/^\s+/, ' ').trim()).filter(line => line);
                const cleaned = cleanedLines.join(' ').trim();
                return `<p${attrs || ''}>${cleaned}</p>`;
            }
            return match;
        });
    }

    htmlToMarkdown(html) {
        // Pre-process to clean up paragraph whitespace
        html = this.cleanParagraphWhitespace(html);
        
        html = html.replace(/<script\b[^<]*(?:(?!<\/script>)<[^<]*)*<\/script>/gis, '');
        html = html.replace(/<style\b[^<]*(?:(?!<\/style>)<[^<]*)*<\/style>/gis, '');
        html = html.replace(/<!--[\s\S]*?-->/g, '');
        html = html.replace(/<nav\b[\s\S]*?<\/nav>/gi, '');
        html = html.replace(/<div[^>]*class=["']manuscript-section[^"']*["'][^>]*data-section=["']([^"']*)["'][^>]*>/gi, '\n\n## $1\n\n');

        html = html.replace(/<pre[^>]*>\s*<code(?: class=["']language-([^"']+)["'])?[^>]*>([\s\S]*?)<\/code>\s*<\/pre>/gi, (match, lang, code) => {
            const language = (lang || '').trim();
            const decodedCode = this.decodeEntities(code).replace(/\n+$/g, '');
            return `\n\n@@@CODEBLOCK_START:${language}@@@\n${decodedCode}\n@@@CODEBLOCK_END@@@\n\n`;
        });

        html = html.replace(/<pre[^>]*>([\s\S]*?)<\/pre>/gi, (match, content) => {
            const decoded = this.decodeEntities(content).replace(/\n+$/g, '');
            return `\n\n@@@CODEBLOCK_START:@@@\n${decoded}\n@@@CODEBLOCK_END@@@\n\n`;
        });

        html = html.replace(/<table[^>]*>[\s\S]*?<\/table>/gi, (match) => this.tableToMarkdown(match));

        // Helper to clean heading content
        const cleanHeading = (content) => content.replace(/^\s+/g, '').replace(/\s+$/g, '');
        
        html = html.replace(/<h1[^>]*>([\s\S]*?)<\/h1>/gi, (m, c) => `\n# ${cleanHeading(c)}\n\n`);
        html = html.replace(/<h2[^>]*>([\s\S]*?)<\/h2>/gi, (m, c) => `\n## ${cleanHeading(c)}\n\n`);
        html = html.replace(/<h3[^>]*>([\s\S]*?)<\/h3>/gi, (m, c) => `\n### ${cleanHeading(c)}\n\n`);
        html = html.replace(/<h4[^>]*>([\s\S]*?)<\/h4>/gi, (m, c) => `\n#### ${cleanHeading(c)}\n\n`);
        html = html.replace(/<h5[^>]*>([\s\S]*?)<\/h5>/gi, (m, c) => `\n##### ${cleanHeading(c)}\n\n`);
        html = html.replace(/<h6[^>]*>([\s\S]*?)<\/h6>/gi, (m, c) => `\n###### ${cleanHeading(c)}\n\n`);

        html = html.replace(/<figcaption[^>]*>([\s\S]*?)<\/figcaption>/gi, '\n\n    $1\n');
        html = html.replace(/<figure[^>]*>([\s\S]*?)<\/figure>/gi, '\n\n$1\n\n');

        html = html.replace(/<img[^>]+>/gi, (tag) => {
            const srcMatch = tag.match(/src=["']([^"']+)["']/i);
            const altMatch = tag.match(/alt=["']([^"']*)["']/i);
            const src = srcMatch ? srcMatch[1] : '';
            const alt = altMatch ? altMatch[1] : '';
            return src ? `\n![${alt}](${src})\n` : '';
        });

        html = html.replace(/<blockquote[^>]*>([\s\S]*?)<\/blockquote>/gi, '\n> $1\n\n');
        html = html.replace(/<p(?:\s[^>]*)?>([\s\S]*?)<\/p>/gi, (match, content) => {
            const trimmed = content.replace(/^\s+/g, '').replace(/\s+$/g, '');
            return trimmed + '\n\n';
        });

        html = html.replace(/<(strong|b)[^>]*>([\s\S]*?)<\/(strong|b)>/gi, '**$2**');
        html = html.replace(/<(em|i)[^>]*>([\s\S]*?)<\/(em|i)>/gi, '*$2*');
        html = html.replace(/<code[^>]*>([\s\S]*?)<\/code>/gi, '`$1`');

        html = html.replace(/<li[^>]*>([\s\S]*?)<\/li>/gi, '- $1\n');
        html = html.replace(/<br\s*\/?>/gi, '\n');
        html = html.replace(/<hr\s*\/?>/gi, '\n---\n');

        html = html.replace(/<\/?[a-zA-Z][^>]*>/g, '');
        html = this.decodeEntities(html);
        html = html.replace(/@@@CODEBLOCK_START:([^@]*)@@@[\r\n]+([\s\S]*?)[\r\n]+@@@CODEBLOCK_END@@@/g, (match, lang, code) => {
            const language = lang.trim();
            return `\n\n\`\`\`${language}\n${code}\n\`\`\`\n\n`;
        });

        // Collapse multiple newlines and trim
        html = html.replace(/\n{3,}/g, '\n\n').trim();
        
        // Final pass: remove leading indentation from all lines (4+ spaces)
        // This handles any residual indentation from HTML structure
        html = html.split('\n').map(line => line.replace(/^\s{4,}/, '')).join('\n');
        
        return html;
    }

    async convertSiteToMarkdown() {
        console.log('🔄 Converting HTML site to markdown (TEP-J0437)...');
        try {
            const manifestPath = path.join(__dirname, 'manifest.json');
            if (!fs.existsSync(manifestPath)) throw new Error('manifest.json not found.');
            const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
            const sections = manifest.sections.sort((a, b) => a.order - b.order);

            // Load citation metadata for header
            const citationPath = path.join(__dirname, '..', 'CITATION.cff');
            let author = 'Matthew Lukin Smawfield';
            let version = 'v0.1 (Sintra)';
            let dateReleased = '2026-05-17';
            let doi = '';
            
            if (fs.existsSync(citationPath)) {
                const citationData = yaml.load(fs.readFileSync(citationPath, 'utf8'));
                if (citationData.authors && citationData.authors[0]) {
                    const firstAuthor = citationData.authors[0];
                    author = `${firstAuthor['given-names']} ${firstAuthor['family-names']}`;
                }
                version = citationData.version || version;
                const rawDate = citationData['date-released'] || dateReleased;
                // Format date as "DD Month YYYY"
                const dateObj = new Date(rawDate);
                const options = { day: 'numeric', month: 'long', year: 'numeric' };
                dateReleased = dateObj.toLocaleDateString('en-GB', options);
                doi = citationData.doi || '';
            }

            let allHtml = '';
            for (const section of sections) {
                const componentPath = path.join(__dirname, 'components', section.file);
                if (fs.existsSync(componentPath)) {
                    let html = fs.readFileSync(componentPath, 'utf8');
                    allHtml += `\n<!-- SECTION: ${section.title} -->\n${html}\n`;
                    console.log(`  ✓ ${section.file} (${(html.length / 1024).toFixed(1)} KB)`);
                } else {
                    console.warn(`  ⚠ Missing: ${section.file}`);
                }
            }

            console.log(`  Total HTML: ${(allHtml.length / 1024).toFixed(1)} KB`);
            const markdown = this.htmlToMarkdown(allHtml);
            
            // Add title and header metadata from citation
            const title = manifest.title || 'Untitled';
            const header = `# ${title}
**${author}**
Version: ${version}
First published: ${dateReleased}${doi ? `\nDOI: ${doi}` : ''}

---

`;
            const finalMarkdown = header + markdown;
            
            // Generate filename: {paper_number}-TEP-J0437-{version}-{codename}.md
            const paperNum = manifest.paper_number || 16;
            const versionMatch = version.match(/(v[\d.]+)\s*\(([^)]+)\)/);
            const versionClean = versionMatch ? versionMatch[1] : version;
            const codename = versionMatch ? versionMatch[2] : 'Sintra';
            const filename = `${paperNum}-TEP-J0437-${versionClean}-${codename}.md`;
            const outputPath = path.join(__dirname, '..', filename);
            
            fs.writeFileSync(outputPath, finalMarkdown, 'utf8');
            console.log(`✅ Markdown saved to: ${outputPath} (${(finalMarkdown.length / 1024).toFixed(1)} KB)`);
        } catch (error) {
            console.error('❌ Markdown conversion failed:', error.message);
        }
    }
}

if (require.main === module) {
    const c = new HTMLToMarkdownConverter();
    c.convertSiteToMarkdown();
}
module.exports = { HTMLToMarkdownConverter };
