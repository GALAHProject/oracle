
from __future__ import print_function

import os
import sys

import github


def get_last_commit_in_pull_request(auth_token, repo_slug, pull_request):
    """
    Get the last commit in a given repository's pull request.
    """

    repo = auth_token.get_repo(repo_slug)
    pr = repo.get_pull(pull_request)

    # Get last commit. This is probably the worst way possible.
    for commit in pr.get_commits():
        continue

    return commit


def get_commit_info(auth_token, repo_slug, pull_request, context):
    """
    Check if this pull request has a pending state with the same context.
    """

    commit = get_last_commit_in_pull_request(auth_token, repo_slug, pull_request)

    context_states = []
    for status in commit.get_statuses():
        if status.raw_data["context"] == context: 
            context_states.append(status.raw_data["state"])

    context_states = list(set(context_states))

    return (commit, context_states)



if __name__ == "__main__":

    # Check for credentials
    if "GH_TOKEN" not in os.environ:
        sys.exit(1)

    # Check that we are on a pull request
    if os.environ["TRAVIS_PULL_REQUEST"] == "false":
        sys.exit(0)

    gh = github.Github(os.environ["GH_TOKEN"])
    context = "science-verification/gaia-benchmarks"
    repo_slug = os.environ["TRAVIS_REPO_SLUG"]
    pull_request = int(os.environ["TRAVIS_PULL_REQUEST"])

    commit, states = get_commit_info(gh, repo_slug, pull_request, context)

    if len(states) == 0:
        # Entry run
        r = commit.create_status("pending", target_url="http://astrowizici.st",
            description="Analysing benchmark stars", context=context)
        sys.exit(0)

    else:
        # Exit run
        r = commit.create_status("success", target_url="http://astrowizici.st",
            description="Gaia benchmark stars within acceptable uncertainties",
            context=context)
        sys.exit(0)


