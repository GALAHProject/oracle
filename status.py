# coding: utf-8

""" Update GitHub pull requests with science verification statuses """ 

from __future__ import print_function

import os
import sys
from collections import deque 

import github


def get_last_commit_in_pull_request(auth_token, repo_slug, pull_request):
    """
    Get the last commit in a given repository's pull request.
    """

    pr = auth_token.get_repo(repo_slug).get_pull(pull_request)

    # Get last commit. This is probably the second worst way possible.
    return deque(pr.get_commits(), maxlen=1).pop()


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
    print("Checking to see if we are on a pull request: {0}".format(
        os.environ.get("TRAVIS_PULL_REQUEST", None)))
    if os.environ.get("TRAVIS_PULL_REQUEST", "false") == "false":
        print("Exiting..")
        sys.exit(0)

    print("We are on pull request #{0}, so we have work to do..".format(
        os.environ["TRAVIS_PULL_REQUEST"]))
    gh = github.Github(os.environ["GH_TOKEN"])
    context = "science-verification/gaia-benchmarks"
    repo_slug = os.environ["TRAVIS_REPO_SLUG"]
    pull_request = int(os.environ["TRAVIS_PULL_REQUEST"])

    commit, states = get_commit_info(gh, repo_slug, pull_request, context)

    if len(states) == 0:
        print("This is the entry run. Set status to pending.")
        # Entry run
        r = commit.create_status("pending", target_url="http://astrowizici.st",
            description="Analysing benchmark stars", context=context)
        sys.exit(0)

    else:
        # Exit run
        print("Current state of PR #{0} is {1}, checking for science results..".format(
            pull_request, "|".join(states)))
        # Was any science actually done?
        results_filename = "science.log"
        if os.path.exists(results_filename):
            print("SCIENCE LOG CONTAINS")
            with open(results_filename, "r") as fp:
                contents = fp.read()

            pr = auth_token.get_repo(repo_slug).get_pull(pull_request)
            new_comment = pr.create_issue_comment("Well:\n" + "".join(contents))

            # [TODO] Parse the log/similar for results and make checks.
            # Then set as either success/failure
            r = commit.create_status("success", target_url="http://astrowizici.st",
                description="Gaia benchmarks within mutual uncertainties",
                context=context)

        else:
            print("No results were found in {0}".format(results_filename))

            pr = auth_token.get_repo(repo_slug).get_pull(pull_request)
            new_comment = pr.create_issue_comment("Could not find `{0}`".format(
                results_filename))

            r = commit.create_status("error", target_url="http://astrowizici.st",
                description="An error occurred and no results were found",
                context=context)

        sys.exit(0)


