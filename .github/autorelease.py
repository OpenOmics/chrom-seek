#!/usr/bin/env python3
import subprocess
import json
import argparse
import sys
import os
import re
from urllib.parse import urlparse
from github import Github, GithubException
from github import Auth

# statics
RELEASE_NOTES_DELIMITER = "## Release notes"
RELEASE_OPTIONS_DELIMITER = "## Release Options"


def create_github_tag(repo_url, token, tag_name, message="", target_sha=None):
    """
    Create a tag on a GitHub repository
    
    Args:
        repo_url: Full GitHub repository URL (e.g., "https://github.com/owner/repo")
        token: GitHub personal access token
        tag_name: Name of the tag (e.g., "v1.0.0")
        message: Optional tag message/annotation
        target_sha: Specific commit SHA to tag (default: latest commit on default branch)
    
    Returns:
        GitTag object
    """
    # Extract owner/repo from URL
    match = re.search(r'github\.com[:/]([^/]+)/([^/.]+)', repo_url)
    if not match:
        raise ValueError(f"Invalid GitHub URL: {repo_url}")
    
    owner = match.group(1)
    repo_name = match.group(2)
    full_repo = f"{owner}/{repo_name}"
    
    # Initialize GitHub client
    g = Github(auth=Auth.Token(args.token))
    repo = g.get_repo(full_repo)
    
    # Get target commit SHA if not provided
    if target_sha is None:
        # Use the latest commit on the default branch
        default_branch = repo.default_branch
        target_sha = repo.get_branch(default_branch).commit.sha
    
    # Get the commit object
    commit = repo.get_commit(target_sha)
    
    # Create a lightweight tag (just a reference)
    if not message:
        ref = repo.create_git_ref(
            ref=f"refs/tags/{tag_name}",
            sha=target_sha
        )
        # print(f"Lightweight tag '{tag_name}' created at {target_sha}")
        return ref
    
    # Create an annotated tag (with message and metadata)
    else:
        # First create the tag object
        tag = repo.create_git_tag(
            tag=tag_name,
            message=message,
            object=target_sha,
            type="commit"
        )
        
        # Then create the reference to the tag
        ref = repo.create_git_ref(
            ref=f"refs/tags/{tag_name}",
            sha=tag.sha
        )
        # print(f"Annotated tag '{tag_name}' created at {target_sha}")
        # print(f"Message: {message}")
        return tag


def create_tag_and_release(repo_url, token, tag_name, release_name, release_body="", tag_message=""):
    """
    Create both a tag and a release in one operation
    
    Args:
        repo_url: Full GitHub repository URL
        token: GitHub personal access token
        tag_name: Name of the tag (e.g., "v1.0.0")
        release_name: Name of the release
        release_body: Release description/notes
        tag_message: Optional tag annotation message
    
    Returns:
        GitRelease object
    """
    # Extract owner/repo from URL
    match = re.search(r'github\.com[:/]([^/]+)/([^/.]+)', repo_url)
    if not match:
        raise ValueError(f"Invalid GitHub URL: {repo_url}")
    
    owner = match.group(1)
    repo_name = match.group(2)
    full_repo = f"{owner}/{repo_name}"
    
    # Initialize GitHub client
    g = Github(auth=Auth.Token(args.token))
    repo = g.get_repo(full_repo)
    
    # Create release (this automatically creates the tag)
    release = repo.create_git_release(
        tag=tag_name,
        name=release_name,
        message=release_body,
        draft=False,
        prerelease=False
    )
    
    # print(f"Tag and release '{tag_name}' created")
    # print(f"Release URL: {release.html_url}")
    return release


def increment_version(version_string, options):
    """Increment major, minor, or patch version"""
    type = None
    assert sum(options.values()) == 1, 'Multiple version types specified'
    parts = list(map(int, version_string.split('.')))

    if options['major']:
        parts[0] += 1
        parts[1] = 0
        parts[2] = 0
    elif options['minor']:
        parts[1] += 1
        parts[2] = 0
    elif options['patch']:
        parts[2] += 1

    return '.'.join(map(str, parts))


def get_latest_version(repo_url, token=None):
    """
    Get the latest release version from a GitHub repo.
    Returns version as a string (e.g., "1.2.3").
    If no releases exist, returns "0.0.0".
    
    Args:
        repo_url: GitHub repo URL (e.g., "https://github.com/owner/repo")
        token: Optional GitHub personal access token (for higher rate limits)
    
    Returns:
        str: Version string
    """
    # Extract owner/repo from URL
    if repo_url.startswith("https://github.com/"):
        repo_path = repo_url.replace("https://github.com/", "").rstrip('/')
        # Remove .git suffix if present
        repo_path = repo_path.replace(".git", "")
    else:
        repo_path = repo_url
    
    # Initialize GitHub client
    # If no token provided, uses unauthenticated requests (lower rate limit)
    g = Github(auth=Auth.Token(args.token))
    
    try:
        # Get the repository
        repo = g.get_repo(repo_path)
        
        # Get the latest release
        latest_release = repo.get_latest_release()
        
        # Get the tag name
        tag = latest_release.tag_name
        
        # Remove 'v' prefix if present (e.g., "v1.2.3" -> "1.2.3")
        version = tag.lstrip('v')
        return version
        
    except GithubException as e:
        # Handle "Not Found" error (no releases exist)
        if e.status == 404:
            return "0.0.0"
        else:
            # Re-raise other exceptions
            raise
    except Exception as e:
        # Handle any other errors
        # print(f"Error getting latest version: {e}")
        return "0.0.0"
    

def get_pr_body(url):
    result = subprocess.run(
        ["gh", "pr", "view", url, "--json", "body"],
        capture_output=True,
        text=True,
        check=True
    )
    body_json = json.loads(result.stdout)["body"]
    return [line for line in body_json.splitlines() if line != '']


def get_notes_and_options(pr_url):
    """
    Parse release options from markdown checklist.
    
    Args:
        url: pull request url
    
    Returns:
        notes: list[str]: lines of release notes
        option: dict[bool]: which type of release or skip
    """
    global RELEASE_NOTES_DELIMITER, RELEASE_OPTIONS_DELIMITER
    content = get_pr_body(pr_url)
    notes_title = content.index(RELEASE_NOTES_DELIMITER)
    options_title = content.index(RELEASE_OPTIONS_DELIMITER)
    notes_content = content[notes_title+1:options_title]
    options_content = content[options_title+1:]

    return "\n".join(notes_content), parse_release_options(options_content)


def parse_release_options(lines):
    """
    Parse release options from markdown checklist.
    
    Args:
        lines: List of strings in markdown checkbox format
    
    Returns:
        dict: Keys are release types, values are booleans
    """
    result = {
        'major': False,
        'minor': False,
        'patch': False,
        'skip': False
    }
    
    for line in lines:
        line_lower = line.lower()
        is_checked = '[x]' in line_lower
        
        if 'major' in line_lower:
            result['major'] = is_checked
        elif 'minor' in line_lower:
            result['minor'] = is_checked
        elif 'patch' in line_lower:
            result['patch'] = is_checked
        elif 'skip' in line_lower:
            result['skip'] = is_checked
    
    return result


def read_version_file(version_file_path="VERSION"):
    """
    Read the version from the VERSION file in the repository root.
    
    Args:
        version_file_path: Path to the VERSION file (default: "VERSION")
    
    Returns:
        str: Version string from the file, stripped of whitespace
    """
    try:
        with open(version_file_path, 'r') as f:
            return f.read().strip()
    except FileNotFoundError:
        raise FileNotFoundError(f"VERSION file not found at {version_file_path}")
    except Exception as e:
        raise Exception(f"Error reading VERSION file: {e}")


def check_version_file(args):
    """
    Check if the VERSION file contains the expected version based on PR options.
    This should be run before PR merge.
    
    Args:
        args: Parsed command-line arguments
    
    Returns:
        bool: True if version is correct, False otherwise
    
    Raises:
        SystemExit: If version check fails
    """
    def extract_repo_url_from_pr(pr_url):
        return pr_url.rstrip('/').split('/pull/')[0]

    PR_URL = args.url
    REPO_URL = extract_repo_url_from_pr(PR_URL)
    
    # Get PR options
    pr_notes, pr_options = get_notes_and_options(PR_URL)
    
    # If skip is selected, version file doesn't need to be updated
    if pr_options['skip']:
        print("✓ Release skipped - VERSION file check not required")
        return True
    
    # Verify exactly one release type is selected
    release_types = [k for k, v in pr_options.items() if k != 'skip' and v]
    if len(release_types) != 1:
        print(f"✗ ERROR: Must select exactly one release type (found: {len(release_types)})")
        sys.exit(1)
    
    # Get current version from GitHub releases
    current_version = get_latest_version(REPO_URL, token=args.token)
    
    # Calculate expected version based on PR options
    expected_version = increment_version(current_version, pr_options)
    
    # Read VERSION file
    try:
        version_file_content = read_version_file(args.version_file)
    except Exception as e:
        print(f"✗ ERROR: {e}")
        sys.exit(1)
    
    # Compare versions
    if version_file_content == expected_version:
        release_type = release_types[0]
        print(f"✓ VERSION file is correct!")
        print(f"  Current version: {current_version}")
        print(f"  Expected version ({release_type}): {expected_version}")
        print(f"  VERSION file: {version_file_content}")
        return True
    else:
        release_type = release_types[0]
        print(f"✗ ERROR: VERSION file does not match expected version!")
        print(f"  Current version: {current_version}")
        print(f"  Expected version ({release_type}): {expected_version}")
        print(f"  VERSION file: {version_file_content}")
        print(f"\nPlease update the VERSION file to: {expected_version}")
        sys.exit(1)


def main(args):
    # If check-version flag is set, run version check and exit
    if args.check_version:
        check_version_file(args)
        return
    
    def extract_repo_url_from_pr(pr_url):
        return pr_url.rstrip('/').split('/pull/')[0]

    PR_URL = args.url
    REPO_URL = extract_repo_url_from_pr(PR_URL)
    # get the full content of the pull request
    pr_data_full = get_pr_body(PR_URL)
    # slice content into notes & options
    pr_notes, pr_options = get_notes_and_options(PR_URL)
    if pr_options['skip']:
        # early exit if no release
        return
    # get current version
    current_version = get_latest_version(REPO_URL)
    # increment current version based on options
    new_version = increment_version(current_version, pr_options)

    create_tag_and_release(
        repo_url=args.url,
        token=args.token,
        tag_name=f"v{new_version}",
        release_name=f"Version {new_version}",
        release_body=pr_notes
    )
    print(f"v{new_version}")
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Automate GitHub releases and version checks for pull requests"
    )
    
    parser.add_argument(
        "url",
        type=str,
        default=os.getenv('PR_URL'),
        help="Valid URL to process (default from PR_URL env var)"
    )

    parser.add_argument(
        "-n", "--pr-number",
        type=str,
        default=os.getenv('PR_NUMBER'),
        help="Pull request number"
    )

    parser.add_argument(
        "-t", "--token",
        type=str,
        default=os.getenv('GITHUB_TOKEN'),
        help="GitHub token"
    )
    
    parser.add_argument(
        "--check-version",
        action="store_true",
        help="Check if VERSION file matches expected version (pre-merge check)"
    )
    
    parser.add_argument(
        "--version-file",
        type=str,
        default="VERSION",
        help="Path to the VERSION file (default: VERSION)"
    )
    
    args = parser.parse_args()
    main(args)
