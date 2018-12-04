library(scatterplot3d)
library(rgl)
library(matlib)
library(plotly)

trials <- c("data/ur5testresult_fullspeed_payload1.6lb_1.csv",
            "data/ur5testresult_fullspeed_payload1.6lb_2.csv",
            "data/ur5testresult_fullspeed_payload1.6lb_3.csv",
            "data/ur5testresult_fullspeed_payload4.5lb_1.csv",
            "data/ur5testresult_fullspeed_payload4.5lb_2.csv",
            "data/ur5testresult_fullspeed_payload4.5lb_3.csv",
            "data/ur5testresult_halfspeed_payload1.6lb_1.csv",
            "data/ur5testresult_halfspeed_payload1.6lb_2.csv",
            "data/ur5testresult_halfspeed_payload1.6lb_3.csv",
            "data/ur5testresult_halfspeed_payload4.5lb_1.csv",
            "data/ur5testresult_halfspeed_payload4.5lb_2.csv",
            "data/ur5testresult_halfspeed_payload4.5lb_3.csv",
            "data/ur5testresult_coldstart_halfspeed_payload4.5lb_1.csv",
            "data/ur5testresult_coldstart_halfspeed_payload4.5lb_2.csv",
            "data/ur5testresult_coldstart_fullspeed_payload4.5lb_3.csv",
            "data/ur5testresult_coldstart_fullspeed_payload4.5lb_1.csv",
            "data/ur5testresult_coldstart_fullspeed_payload4.5lb_2.csv",
            "data/ur5testresult_coldstart_fullspeed_payload4.5lb_3.csv")

df <- data.frame(time = seq(0, 47.992, length = 6000))

for (name in trials) {
  ur5 <- read.csv(file = name, header = FALSE)
  ur5.headers <- scan("data/ur5_header.csv", sep = ',', what = "", quiet = TRUE)
  colnames(ur5) <- ur5.headers
  x <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (x)`
  y <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (y)`
  z <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (z)`
  rx <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (rx)`
  ry <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (ry)`
  rz <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (rz)`

  x <- gsub("[^0-9\\.]", "", x) # cleans non-numeric characters
  x <- as.numeric(x)

  rz <- gsub("[^0-9\\.]", "", rz) # cleans non-numeric characters
  rz <- as.numeric(rz)

  rx <- rx[1:6000]
  ry <- ry[1:6000]
  rz <- rz[1:6000]

  df$placeholder_name = rx
  names(df)[names(df) == "placeholder_name"] <- paste(name, "rx")

  df$placeholder_name = ry
  names(df)[names(df) == "placeholder_name"] <- paste(name, "ry")

  df$placeholder_name = rz
  names(df)[names(df) == "placeholder_name"] <- paste(name, "rz")

}

meanrx1.6 <- (df$`data/ur5testresult_fullspeed_payload1.6lb_1.csv rx` +
              df$`data/ur5testresult_fullspeed_payload1.6lb_2.csv rx` +
              df$`data/ur5testresult_fullspeed_payload1.6lb_3.csv rx`)/3

meanry1.6 <- (df$`data/ur5testresult_fullspeed_payload1.6lb_1.csv ry` +
              df$`data/ur5testresult_fullspeed_payload1.6lb_2.csv ry` +
              df$`data/ur5testresult_fullspeed_payload1.6lb_3.csv ry`)/3

meanrz1.6 <- (df$`data/ur5testresult_fullspeed_payload1.6lb_1.csv rz` +
              df$`data/ur5testresult_fullspeed_payload1.6lb_2.csv rz` +
              df$`data/ur5testresult_fullspeed_payload1.6lb_3.csv rz`)/3

meanrx4.5 <- (df$`data/ur5testresult_fullspeed_payload4.5lb_1.csv rx` +
              df$`data/ur5testresult_fullspeed_payload4.5lb_2.csv rx` +
              df$`data/ur5testresult_fullspeed_payload4.5lb_3.csv rx`)/3

meanry4.5 <- (df$`data/ur5testresult_fullspeed_payload4.5lb_1.csv ry` +
              df$`data/ur5testresult_fullspeed_payload4.5lb_2.csv ry` +
              df$`data/ur5testresult_fullspeed_payload4.5lb_3.csv ry`)/3

meanrz4.5 <- (df$`data/ur5testresult_fullspeed_payload4.5lb_1.csv rz` +
              df$`data/ur5testresult_fullspeed_payload4.5lb_2.csv rz` +
              df$`data/ur5testresult_fullspeed_payload4.5lb_3.csv rz`)/3

mean.sq.error.x <- mean((meanrx1.6 - meanrx4.5) ^ 2)
mean.sq.error.y <- mean((meanry1.6 - meanry4.5) ^ 2)
mean.sq.error.z <- mean((meanrz1.6 - meanrz4.5) ^ 2)


scatterplot3d(x, y, z)

plot_ly(
  type = "cone",
  x = c(1, 0, 0),
  y = c(0, 1, 0),
  z = c(0, 0, 1),
  u = c(-mean.sq.error.x, 0, 0),
  v = c(0, -mean.sq.error.y, 0),
  w = c(0, 0, -mean.sq.error.z),
  sizemode = "absolute",
  cmin = 0,
  cmax = 0.03,


  anchor = "base")%>%
  add_trace(type = 'scatter3d', x = c(0,1), y = c(0,0), z = c(0,0), mode = 'lines',
            opacity = 1, line = list(width = 6, reverscale = FALSE), showlegend = FALSE)%>%
  add_trace(type = 'scatter3d', x = c(0,0), y = c(0,1), z = c(0,0), mode = 'lines',
            opacity = 1, line = list(width = 6, reverscale = FALSE), showlegend = FALSE)%>%
  add_trace(type = 'scatter3d', x = c(0,0), y = c(0,0), z = c(0,1), mode = 'lines',
            opacity = 1, line = list(width = 6, reverscale = FALSE), showlegend = FALSE)%>%
  layout(
    title = "Mean Squared Error in Rotation Vector: 4.5lb vs 1.6lb load",
    scene = list(
      xaxis = list(
        title = "Y",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      yaxis = list(
        title = "X",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      zaxis = list(
        title = "Z",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      )
    ))


meanrxcold <- (df$`data/ur5testresult_coldstart_fullspeed_payload4.5lb_1.csv rx` +
               df$`data/ur5testresult_coldstart_fullspeed_payload4.5lb_2.csv rx` +
               df$`data/ur5testresult_coldstart_fullspeed_payload4.5lb_3.csv rx`)/3

meanrycold <- (df$`data/ur5testresult_coldstart_fullspeed_payload4.5lb_1.csv ry` +
               df$`data/ur5testresult_coldstart_fullspeed_payload4.5lb_2.csv ry` +
               df$`data/ur5testresult_coldstart_fullspeed_payload4.5lb_3.csv ry`)/3

meanrzcold <- (df$`data/ur5testresult_coldstart_fullspeed_payload4.5lb_1.csv rz` +
               df$`data/ur5testresult_coldstart_fullspeed_payload4.5lb_2.csv rz` +
               df$`data/ur5testresult_coldstart_fullspeed_payload4.5lb_3.csv rz`)/3

meanrxwarm <- (df$`data/ur5testresult_fullspeed_payload4.5lb_1.csv rx` +
              df$`data/ur5testresult_fullspeed_payload4.5lb_2.csv rx` +
              df$`data/ur5testresult_fullspeed_payload4.5lb_3.csv rx`)/3

meanrywarm <- (df$`data/ur5testresult_fullspeed_payload4.5lb_1.csv ry` +
              df$`data/ur5testresult_fullspeed_payload4.5lb_2.csv ry` +
              df$`data/ur5testresult_fullspeed_payload4.5lb_3.csv ry`)/3

meanrzwarm <- (df$`data/ur5testresult_fullspeed_payload4.5lb_1.csv rz` +
              df$`data/ur5testresult_fullspeed_payload4.5lb_2.csv rz` +
              df$`data/ur5testresult_fullspeed_payload4.5lb_3.csv rz`)/3

mean.sq.error.x <- mean((meanrxcold - meanrxwarm) ^ 2)
mean.sq.error.y <- mean((meanrycold - meanrywarm) ^ 2)
mean.sq.error.z <- mean((meanrzcold - meanrzwarm) ^ 2)

plot_ly(
  type = "cone",
  x = c(1, 0, 0),
  y = c(0, 1, 0),
  z = c(0, 0, 1),
  u = c(-mean.sq.error.x, 0, 0),
  v = c(0, -mean.sq.error.y, 0),
  w = c(0, 0, -mean.sq.error.z),
  sizemode = "absolute",
  cmin = 0,
  cmax = 0.03,


  anchor = "base")%>%
  add_trace(type = 'scatter3d', x = c(0,1), y = c(0,0), z = c(0,0), mode = 'lines',
            opacity = 1, line = list(width = 6, reverscale = FALSE), showlegend = FALSE)%>%
  add_trace(type = 'scatter3d', x = c(0,0), y = c(0,1), z = c(0,0), mode = 'lines',
            opacity = 1, line = list(width = 6, reverscale = FALSE), showlegend = FALSE)%>%
  add_trace(type = 'scatter3d', x = c(0,0), y = c(0,0), z = c(0,1), mode = 'lines',
            opacity = 1, line = list(width = 6, reverscale = FALSE), showlegend = FALSE)%>%
  layout(
    title = "Mean Squared Error in Rotation Vector: Warm vs Cold Start",
    scene = list(
      xaxis = list(
        title = "Y",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      yaxis = list(
        title = "X",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      zaxis = list(
        title = "Z",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      )
    ))
