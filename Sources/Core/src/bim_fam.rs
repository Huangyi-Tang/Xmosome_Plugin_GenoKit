//! Parsers and serializers for PLINK .fam (samples) and .bim (variants) files.

use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FamRecord {
    pub fid: String,
    pub iid: String,
    pub pid: String,
    pub mid: String,
    pub sex: String,
    pub pheno: String,
}

impl FamRecord {
    pub fn parse_line(line: &str) -> Option<Self> {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            return None;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 6 {
            return None;
        }
        Some(Self {
            fid: parts[0].to_string(),
            iid: parts[1].to_string(),
            pid: parts[2].to_string(),
            mid: parts[3].to_string(),
            sex: parts[4].to_string(),
            pheno: parts[5].to_string(),
        })
    }

    pub fn to_fam_line(&self) -> String {
        format!("{} {} {} {} {} {}", self.fid, self.iid, self.pid, self.mid, self.sex, self.pheno)
    }
}

impl fmt::Display for FamRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_fam_line())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BimRecord {
    pub chr: String,
    pub id: String,
    pub cm: String,
    pub pos: u64,
    pub a1: String,
    pub a2: String,
}

impl BimRecord {
    pub fn parse_line(line: &str) -> Option<Self> {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            return None;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() < 6 {
            return None;
        }
        let pos = parts[3].parse::<u64>().unwrap_or(0);
        Some(Self {
            chr: parts[0].to_string(),
            id: parts[1].to_string(),
            cm: parts[2].to_string(),
            pos,
            a1: parts[4].to_string(),
            a2: parts[5].to_string(),
        })
    }

    pub fn to_bim_line(&self) -> String {
        format!("{}\t{}\t{}\t{}\t{}\t{}", self.chr, self.id, self.cm, self.pos, self.a1, self.a2)
    }
}

impl fmt::Display for BimRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_bim_line())
    }
}

pub fn parse_fam(content: &str) -> Vec<FamRecord> {
    content
        .lines()
        .filter_map(FamRecord::parse_line)
        .collect()
}

pub fn parse_bim(content: &str) -> Vec<BimRecord> {
    content
        .lines()
        .filter_map(BimRecord::parse_line)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fam_parsing_and_formatting() {
        let content = "FAM01 IND01 0 0 1 1\nFAM01 IND02 0 0 2 2\n# Comment\n  \n";
        let fam = parse_fam(content);
        assert_eq!(fam.len(), 2);
        assert_eq!(fam[0].fid, "FAM01");
        assert_eq!(fam[0].iid, "IND01");
        assert_eq!(fam[1].sex, "2");
        assert_eq!(fam[0].to_fam_line(), "FAM01 IND01 0 0 1 1");
    }

    #[test]
    fn test_bim_parsing_and_formatting() {
        let content = "1 rs1001 0 12345 A G\n2 rs1002 0.5 67890 C T\n";
        let bim = parse_bim(content);
        assert_eq!(bim.len(), 2);
        assert_eq!(bim[0].chr, "1");
        assert_eq!(bim[0].id, "rs1001");
        assert_eq!(bim[0].pos, 12345);
        assert_eq!(bim[0].a1, "A");
        assert_eq!(bim[0].a2, "G");
        assert_eq!(bim[0].to_bim_line(), "1\trs1001\t0\t12345\tA\tG");
    }
}
