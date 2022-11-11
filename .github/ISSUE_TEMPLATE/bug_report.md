---
name: Bug report
about: File a bug report
title: ''
labels: ["bug", "triage"]
assignees: ''

---
body:
  - type: markdown
    attributes:
      value: |
        Thanks for using Korg and for taking the time to fill out this bug report!
  - type: textarea
    id: what-happened
    attributes:
      label: What happened?
      description: Also tell us, what did you expect to happen?
      placeholder: Tell us what you see!
      value: 
    validations:
      required: false
  - type: markdown
    attributes:
        value: |
            _You can submit your report without answering the following questions about version numbers.`
  - type: textarea
    id: julia_version
    attributes:
      label: Which version of Julia are you running?
      description: You can check this by entering `julia -v` at your command line. We always recommend using the [current stable release](https://julialang.org/downloads/) of Julia, but please continue to submit your report because we want to make sure Korg works with all versions of Julia.
      value:
    validations:
      required: false
  - type: textarea
    id: korg_version
    attributes:
      label: Which version of Korg are you running?
      description: You can check this by entering `julia -e 'using Pkg; Pkg.status("Korg")'` at a command line.
      value:
  - type: textarea
    id: logs
    attributes:
      label: Relevant log output
      description: Please copy and paste any relevant Julia output. This will be automatically formatted into code, so no need for backticks. If you've already entered the relevant output in the 'What happened?' input then you don't have to replicate it here. 
      render: julia