using Documenter
using jtsk2etrs

makedocs(
    sitename = "jtsk2etrs",
    format = Documenter.HTML(),
    modules = [jtsk2etrs]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
