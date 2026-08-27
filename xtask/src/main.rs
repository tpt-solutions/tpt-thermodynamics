//! Workspace task runner for `tpt-thermodynamics`.
//!
//! Invoked via the `cargo xtask` alias (see `.cargo/config.toml`):
//!
//! ```text
//! cargo xtask fmt              # rustfmt over the whole workspace
//! cargo xtask clippy           # clippy with -D warnings, all features
//! cargo xtask test             # unit/integration tests + doctests
//! cargo xtask deny             # cargo-deny (advisories/bans/licenses/sources)
//! cargo xtask wasm             # cross-check (non-xtask crates) for wasm32
//! cargo xtask check            # fast compile check, all features
//! cargo xtask new-crate <name> # scaffold a new tpt-thermo-* crate
//! cargo xtask all              # everything above, in order
//! ```

use std::fs;
use std::process::{Command, ExitCode};

const WASM_TARGET: &str = "wasm32-unknown-unknown";

const TASKS: &[(&str, &str)] = &[
    ("fmt", "rustfmt over the whole workspace"),
    (
        "clippy",
        "clippy --all-targets --all-features with -D warnings",
    ),
    ("test", "tests + doctests (--all-features)"),
    ("deny", "cargo-deny advisories/bans/licenses/sources"),
    (
        "wasm",
        "cross-check non-xtask crates for wasm32-unknown-unknown",
    ),
    ("check", "fast compile check (--all-features)"),
    (
        "new-crate",
        "scaffold a new tpt-thermo-* crate: new-crate <name>",
    ),
    ("all", "run every task above in order"),
];

fn main() -> ExitCode {
    let mut args = std::env::args().skip(1);
    let task = args.next().unwrap_or_else(|| "help".to_string());
    let ok = match task.as_str() {
        "new-crate" => match args.next() {
            Some(name) => new_crate(&name),
            None => {
                eprintln!("usage: cargo xtask new-crate <name>   (e.g. tpt-thermo-core)");
                false
            }
        },
        other => run(other),
    };
    if ok {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

fn run(task: &str) -> bool {
    match task {
        "help" => {
            print_help();
            true
        }
        "fmt" => cargo(&["fmt", "--all"]),
        "clippy" => cargo(&[
            "clippy",
            "--workspace",
            "--all-targets",
            "--all-features",
            "--",
            "-D",
            "warnings",
        ]),
        "test" => {
            cargo(&["test", "--workspace", "--all-features"])
                && cargo(&["test", "--workspace", "--all-features", "--doc"])
        }
        "deny" => cargo(&["deny", "check"]),
        "wasm" => {
            // `xtask` is a host-only binary (clap) and is excluded; the
            // remaining library crates are checked for the wasm32 target so
            // that the no_std-capable crates (e.g. tpt-thermo-core) stay
            // wasm-buildable.
            cargo(&[
                "check",
                "--workspace",
                "--exclude",
                "xtask",
                "--target",
                WASM_TARGET,
                "--all-features",
            ])
        }
        "check" => cargo(&["check", "--workspace", "--all-features"]),
        "all" => ["fmt", "clippy", "test", "deny", "wasm", "check"]
            .iter()
            .all(|t| run(t)),
        other => {
            eprintln!("unknown task: {other}");
            print_help();
            false
        }
    }
}

/// Run `cargo <args>` in the workspace root, inheriting stdio.
fn cargo(args: &[&str]) -> bool {
    println!("+ cargo {}", args.join(" "));
    Command::new("cargo")
        .args(args)
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn print_help() {
    println!(
        "usage: cargo xtask <task>

tasks:"
    );
    for (name, desc) in TASKS {
        println!("  {name:<10} {desc}");
    }
}

// ---------------------------------------------------------------------------
// Crate scaffolding for the per-phase build-out. Creates
// `crates/<name>/` with a minimal but workspace-lint-clean `Cargo.toml` and
// `src/lib.rs` stub, and registers it in the workspace `members` list.
// ---------------------------------------------------------------------------

fn new_crate(name: &str) -> bool {
    let crate_dir = format!("crates/{name}");
    if std::path::Path::new(&crate_dir).exists() {
        eprintln!("error: {crate_dir} already exists");
        return false;
    }
    let src_dir = format!("{crate_dir}/src");
    if !fs::create_dir_all(&src_dir)
        .map_err(|e| eprintln!("error: {e}"))
        .is_ok()
    {
        return false;
    }

    let cargo_toml = format!(
        "[package]\n\
         name = \"{name}\"\n\
         version = \"0.1.0\"\n\
         edition.workspace = true\n\
         rust-version.workspace = true\n\
         license.workspace = true\n\
         \n\
         [dependencies]\n\
         \n\
         [lints]\n\
         workspace = true\n"
    );
    let lib_rs = format!(
        "//! `{name}` — scaffolded by `cargo xtask new-crate`.\n\
         \n\
         /// Placeholder so the empty crate has at least one public item.\n\
         pub fn placeholder() -> u32 {{\n    \
            0\n\
         }}\n\
         \n\
         #[cfg(test)]\n\
         mod tests {{\n    \
            #[test]\n    \
            fn it_works() {{\n        \
                assert_eq!(placeholder(), 0);\n    \
            }}\n\
         }}\n"
    );

    if let Err(e) = fs::write(format!("{crate_dir}/Cargo.toml"), cargo_toml) {
        eprintln!("error: {e}");
        return false;
    }
    if let Err(e) = fs::write(format!("{src_dir}/lib.rs"), lib_rs) {
        eprintln!("error: {e}");
        return false;
    }

    // Register the member in the workspace manifest if not already present.
    if let Err(e) = register_member(name) {
        eprintln!("warning: could not auto-register member: {e}");
    }

    println!("created {crate_dir}");
    true
}

fn register_member(name: &str) -> Result<(), String> {
    let manifest = "Cargo.toml";
    let content = fs::read_to_string(manifest).map_err(|e| e.to_string())?;
    if content.contains(&format!("\"{name}\"")) || content.contains(&format!("\"crates/{name}\"")) {
        return Ok(());
    }
    // Insert the new member into the `members` array. We key off the existing
    // `members = [` block so the edit is robust to formatting changes.
    let marker = "members = [";
    let idx = content
        .find(marker)
        .ok_or("members array not found in Cargo.toml")?;
    let insert_at = content[idx..]
        .find(']')
        .map(|p| idx + p)
        .ok_or("unterminated members array")?;
    let mut new_content = String::with_capacity(content.len() + 32);
    new_content.push_str(&content[..insert_at]);
    new_content.push_str(&format!("\n    \"crates/{name}\",\n"));
    new_content.push_str(&content[insert_at..]);
    fs::write(manifest, new_content).map_err(|e| e.to_string())?;
    Ok(())
}
