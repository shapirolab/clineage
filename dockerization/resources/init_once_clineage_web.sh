#!/bin/bash
set -e

if [ -f /init_done ]; then
	exit
fi

. $(which activate) py36

cd /clineage
yes "yes" | ./manage.py migrate
cat <<EOF | python manage.py shell
from django.contrib.auth import get_user_model

User = get_user_model()  # get the currently active user model,

User.objects.filter(username='root').exists() or \
    User.objects.create_superuser('root', 'root@example.com', 'root')
EOF

touch /init_done
