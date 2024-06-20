import telegram
import smtplib
from email.mime.text import MIMEText

def SendEmail(From,To,Subject,Content):

  msg = MIMEText(Content)

  msg['Subject'] = Subject
  msg['From'] = From
  msg['To'] = To

  s = smtplib.SMTP('localhost')
  s.sendmail(From, [To], msg.as_string())
  s.quit()

def SendEmailbyGMail(From,To,Subject,Content):

  msg = MIMEText(Content)

  msg['Subject'] = Subject
  msg['From'] = From
  msg['To'] = To

  s = smtplib.SMTP_SSL('smtp.gmail.com', 465)
  s.login('snucms.knu.job.monitor@gmail.com', 'ujuzpiulhdiucnat')
  s.sendmail(From, [To], msg.as_string())
  s.quit()


async def SendTelegramMsg(telegram_id, telegram_token, msg):
    bot = telegram.Bot(telegram_token)
    await bot.send_message(text=text, chat_id=telegram_id)
