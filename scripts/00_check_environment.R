rm(list = ls())

if (!file.exists("DESCRIPTION") || !dir.exists("R") || !dir.exists("scripts")) {
    stop("Run this script from the repository root.")
}

versions <- utils::read.csv(
    "environment/package-versions.csv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)
package_rows <- versions$package != "R"
packages <- versions$package[package_rows]
recommended <- setNames(versions$version[package_rows], packages)
installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)

status <- data.frame(
    package = packages,
    recommended = unname(recommended),
    installed = ifelse(
        installed,
        vapply(packages, function(x) as.character(utils::packageVersion(x)), character(1)),
        NA_character_
    ),
    stringsAsFactors = FALSE
)
print(status, row.names = FALSE)

if (any(!installed)) {
    stop(
        "Missing packages: ",
        paste(packages[!installed], collapse = ", "),
        ". See README.md for installation instructions."
    )
}

recommended_r <- versions$version[versions$package == "R"]
message("R version: ", getRversion(), " (reference run: ", recommended_r, ")")
if (any(status$installed != status$recommended)) {
    message("Some package versions differ from the reference environment; see the table above.")
}
message("Environment check passed.")

