#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Update GitHub pull requests with the science verification status. """ 

from __future__ import print_function

import os
import sys
from collections import deque 

# Note, this requires the forked PyGithub version from 
# https://github.com/andycasey/PyGithub because the standard PyGithub does not
# allow you to provide context when creating a status.
# See https://github.com/jacquev6/PyGithub/pull/289
import github


def get_last_commit_in_pull_request(auth_token, repo_slug, pull_request):
    """ Get the last commit in a given repository's pull request. """

    pr = auth_token.get_repo(repo_slug).get_pull(pull_request)
    # Get last commit. This is probably the second worst way possible.
    return deque(pr.get_commits(), maxlen=1).pop()


def get_commit_info(auth_token, repo_slug, pull_request, context):
    """ Check if the pull request has 'pending' state in the given context. """

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
        print("No GH_TOKEN found in the environment.")
        sys.exit(1)

    # Check that we are on a pull request
    pull_request = os.environ.get("TRAVIS_PULL_REQUEST", None)
    print("Checking if we are on a pull request: {0}".format(pull_request))    
    if pull_request in (None, "false"):
        print("Exiting..")
        sys.exit(0)

    # Create a GH instance and get the information we need.
    gh = github.Github(os.environ["GH_TOKEN"])
    context = "science-verification/gaia-benchmarks"
    repo_slug = os.environ["TRAVIS_REPO_SLUG"]
    pull_request = int(pull_request)

    commit, states = get_commit_info(gh, repo_slug, pull_request, context)

    print("We are on pull request #{0} and the last commit was {1} by @{2} ({3})"\
        .format(pull_request, commit.sha, commit.author.login, commit.author.name))

    if len(states) == 0:

        # Entry run
        print("This is the entry run. Seting status of {0} to 'pending'".format(
            context))
        r = commit.create_status("pending", target_url="http://astrowizici.st",
            description="Analysing benchmark stars", context=context)

    else:
        # Exit run
        print("Current state of PR #{0} is {1}, checking for science results.."\
            .format(pull_request, "|".join(states)))

        # Was any science actually done?
        log_filename = "science/science.log"
        results_filename = "science/results.md"
        if os.path.exists(results_filename):
    
            with open(results_filename, "r") as fp:
                results = fp.read()

            # TODO: Make more things formattable
            formattable = {
                "commit_sha": commit.sha[:10]
            }
            results = results.format(**formattable)

            # Make a comment on the pull request that contains the information
            # from the results file.
            pr = gh.get_repo(repo_slug).get_pull(pull_request)
            new_comment = pr.create_issue_comment(results)

            # [TODO] Parse the log/similar for results and make checks.
            # Then set as either success/failure

            # Update the science-verification state
            r = commit.create_status("success", target_url="http://astrowizici.st",
                description="Gaia benchmarks within mutual uncertainties",
                context=context)

        else:
            print("No results were found in {0}".format(results_filename))

            # Get the log output if it exists.
            if os.path.exists(log_filename):
                with open(log_filename, "r") as fp:
                    log = fp.read()

                else:
                    log = "No log could be found at `{}`".format(log_filename)

            # Make a comment on the pull request that highlights the author of
            # this commit.
            pr = gh.get_repo(repo_slug).get_pull(pull_request)
            new_comment = pr.create_issue_comment("Could not find the science "
                "verification results file `{0}` after analysing with code "
                "commit [`{1}`]({2}) by @{3}:\n\n{4}".format(
                results_filename, commit.sha[:10], commit.url,
                commit.author.login, log))

            # Update the science-verification state
            r = commit.create_status("error", target_url="http://astrowizici.st",
                description="An error occurred and no results were found",
                context=context)

    

